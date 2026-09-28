from pdbfixer import PDBFixer
from openmm.app import *
from openmm import *
from openmm.unit import *
import sys
import time
import argparse


parser_arg = argparse.ArgumentParser()
parser_arg.add_argument(
    "--aa_structure",
    type=str,
    required=True,
    help="aa structure to deform.",
)
parser_arg.add_argument(
    "--output_prefix",
    type=str,
    required=True,
    default="aa",
    help="prefix for the pdb output files after relaxation.",
)


platform = Platform.getPlatformByName("CUDA")

forcefield = ForceField(
    "amber14-all.xml",
    "implicit/gbn2.xml") #Amber ff14SB and implicit solvent


#integrator = LangevinIntegrator(
#    300*kelvin,
#    1/picosecond,
#    0.002*picoseconds)


def aa_move_to_ca(target_pdb_file, output_prefix):
# Load structures
    start_pdb = PDBFile(aa_structure)
    target_pdb = PDBFile(target_pdb_file)

    system = forcefield.createSystem(
    start_pdb.topology,
    nonbondedMethod=NoCutoff,
    constraints=HBonds
    )

# Custom restraint with interpolation parameter lambda
    force = CustomExternalForce(
    "0.5*k*((x-(x0+lambda*(x1-x0)))^2 + (y-(y0+lambda*(y1-y0)))^2 + (z-(z0+lambda*(z1-z0)))^2)"
    )

    force.addGlobalParameter("k", 5000*kilojoule_per_mole/nanometer**2)
    force.addGlobalParameter("lambda", 0.0)

    force.addPerParticleParameter("x0")
    force.addPerParticleParameter("y0")
    force.addPerParticleParameter("z0")

    force.addPerParticleParameter("x1")
    force.addPerParticleParameter("y1")
    force.addPerParticleParameter("z1")

# Store target CA coordinates
    target_ca = {}
    for atom in target_pdb.topology.atoms():
        if atom.name == "CA":
            target_ca[atom.residue.index] = target_pdb.positions[atom.index]

# Add CA restraints
    for atom in start_pdb.topology.atoms():
        if atom.name == "CA":
            res = atom.residue.index
            if res in target_ca:

                start_pos = start_pdb.positions[atom.index].value_in_unit(nanometer)
                target_pos = target_ca[res].value_in_unit(nanometer)

                params = list(start_pos) + list(target_pos)

                force.addParticle(atom.index, params)

    system.addForce(force)

    integrator = LangevinIntegrator(
        300*kelvin,
        1/picosecond,
        0.002*picoseconds)

    simulation = Simulation(start_pdb.topology, system, integrator)
    simulation.context.setPositions(start_pdb.positions)

# Minimize first
    simulation.minimizeEnergy()

# Gradually increase lambda
    n_steps = 20000
    n_windows = 50

    for i in range(n_windows):

        lam = i/(n_windows-1)
        simulation.context.setParameter("lambda", lam)

        simulation.step(n_steps//n_windows)

# Save final structure
    state = simulation.context.getState(getPositions=True)

    with open(f"{output_prefix}_{target_pdb_file}","w") as f:
        PDBFile.writeFile(simulation.topology, state.getPositions(), f)

if __name__ == '__main__':
    #s = time.time()
    args = parser_arg.parse_args()
    aa_structure = args.aa_structure
    output_prefix = args.output_prefix
    #for input_pdb in inputs:
    #    minimize(input_pdb, integrator, platform, output_prefix)
    #print((time.time()-s))

    target_pdbs = [f"cluster_all.c{i}.pdb" for i in range(100)]
    for i in target_pdbs:
        aa_move_to_ca(i, output_prefix)
