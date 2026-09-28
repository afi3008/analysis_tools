from pdbfixer import PDBFixer
from openmm.app import *
from openmm import *
from openmm.unit import *
import sys
import time
import argparse


parser_arg = argparse.ArgumentParser()
#parser_arg.add_argument(
#    "--path_structures",
#    type=str,
#    required=True,
#    help="path to the folder containing the cluster files.",
#)
parser_arg.add_argument(
    "--output_prefix",
    type=str,
    required=True,
    default="relaxed",
    help="prefix for the pdb output files after relaxation.",
)

inputs = [f"aa_cluster_all.c{i}.pdb" for i in range(100)]

platform = Platform.getPlatformByName("CUDA")

forcefield = ForceField(
    "amber14-all.xml",
    "implicit/gbn2.xml") #Amber ff14SB and implicit solvent


def fix_pdb(input_pdb):
    """
    Adds missing atoms to the input pdb-file.
    input_pdb: Input pdb-file.
    return: returns fixed pdb.
    """
    fixer = PDBFixer(filename=input_pdb)
    fixer.findMissingResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()
    fixer.addMissingHydrogens(pH=7.0)
    return fixer


def define_system(top):
    """
    Sets up the system for the minimization/relaxation.
    top: Topology of the pdb-file.
    return: System.
    """
    system = forcefield.createSystem(
        top,
        nonbondedMethod=NoCutoff,
        constraints=None)
    return system


def add_restraints(pdb):
    """
    Add restraints to the heavy atoms of the system. (harmonic potential)
    pdb: fixed pdb-file.
    return: restrained system.
    """
    system = define_system(pdb.topology)
    restraint = CustomExternalForce(
        'k*((x-x0)^2 + (y-y0)^2 + (z-z0)^2)')
    restraint.addGlobalParameter("k", 2000.0 * kilojoule_per_mole / nanometer**2)
    restraint.addPerParticleParameter("x0")
    restraint.addPerParticleParameter("y0")
    restraint.addPerParticleParameter("z0")

    for i, atom in enumerate(pdb.topology.atoms()):
        #print("atom", atom)
        if atom.name == "CA":
            #print("atom_H", atom)
            pos = pdb.positions[i]
            restraint.addParticle(i, pdb.positions[i])

    system.addForce(restraint)
    return system


def minimize_structure(pdb, platform):
    """
    Minimizes/relaxes the sidechains of the Structure. 
    pdb: fixed pdb-file.
    integrator: Langevin integrator.
    platform: Platform to calculate on (CUDA, OpenCL, CPU...).
    return: Minimization Simulation.
    """
    system = add_restraints(pdb)

    integrator = LangevinIntegrator(
    300*kelvin,
    1/picosecond,
    0.002*picoseconds)
    simulation = Simulation(pdb.topology, system, integrator, platform)
    simulation.context.setPositions(pdb.positions)

    # Energy minimization
    simulation.minimizeEnergy(maxIterations=0)  # 0 = until convergence
    return simulation


def get_relaxed(pdb, platform):
    """
    Get new positions of the last step of the minimization simulation.
    return: new positions.
    """
    # Get minimized structure
    simulation = minimize_structure(pdb, platform)
    state = simulation.context.getState(getPositions=True)
    positions = state.getPositions()
    return positions


def write_relax(output, top, positions):
    """
    write out relaxed strucutres.
    return: None.
    """
    with open(output, "w") as f:
        PDBFile.writeFile(top, positions, f)
    
    print(f"Relaxation complete for {output}.")


def minimize(input_pdb, platform, output_prefix):
    """
    Load a pdb, fix the missing atoms and minimize the structure with a restrained on the heavy atoms and save a new structure.
    output_pdb: how is the output called.
    """
    pdb = fix_pdb(input_pdb)
    positions_relax = get_relaxed(pdb, platform)
    write_relax(f"{output_prefix}_{input_pdb}", pdb.topology, positions_relax)


if __name__ == '__main__':
    s = time.time()
    args = parser_arg.parse_args()
    output_prefix = args.output_prefix
    for input_pdb in inputs:
        minimize(input_pdb, platform, output_prefix)
    print((time.time()-s))
