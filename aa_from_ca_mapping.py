from pdbfixer import PDBFixer
from openmm.app import *
from openmm import *
from openmm.unit import *
import sys
import time
import argparse
import numpy as np
import pandas as pd

parser_arg = argparse.ArgumentParser()
#parser_arg.add_argument(
#    "--aa_structures",
#    nargs="+",
#    required=True,
#    help="list of aa structures to deform (bigger cluster aa structures).Increasing numbering! E.g. cluster_4.c0.pdb cluster_4.c1.pdb cluster_4.c2.pdb cluster_4.c3.pdb",
#)
parser_arg.add_argument(
    "--mapping_smaller",
    type=str,
    help="mapping to smaller/more clusters",)
parser_arg.add_argument(
    "--mapping_bigger",
    type=str,
    help="mapping to bigger/less clusters",)
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
def read_mapping(mapping_path):
    """
    Read mapping from clustering (which particle belongs to which cluster).
    """
    cluster_mapping = np.genfromtxt(mapping_path, delimiter=None, dtype=int)
    mapping = cluster_mapping[:, 1]
    return mapping

def map_to_map(map_small, map_big):
    """
    Combine smaller clustering with bigger clustering to dataframe.
    """
    maps = pd.DataFrame({'map_small': map_small, 'map_big': map_big}, columns=['map_small', 'map_big'])
    summary = pd.DataFrame()
    for i in range(maps.map_small.max()+1):
        filtered = maps[maps.map_small == i]
        total = filtered.shape[0]
        #print("total", total)
        list_per_small_cluster = []
        for j in range(maps.map_big.max()+1):
            list_per_small_cluster.append(filtered[filtered.map_big == j].shape[0] / total)
        per_cluster = pd.DataFrame({i: list_per_small_cluster})
        #print("per_cluster",per_cluster)
        summary = pd.concat([summary, per_cluster.T], ignore_index = True)
    return summary

def aa_move_to_ca(aa_structure, target_pdb_file, output_prefix):
    """
    transform aa structure (bigger cluster) to respective CA positions of smaller cluster center.
    """
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
        if aa_struct_idx > -1:
            if atom.name == "CA" and atom.residue.index > 9:
                target_ca[atom.residue.index] = target_pdb.positions[atom.index]
        else:
            if atom.name == "CA":
                target_ca[atom.residue.index] = target_pdb.positions[atom.index]

# Add CA restraints
    for i, atom in enumerate(start_pdb.topology.atoms()):
        if atom.name == "CA":
            res = atom.residue.index
            if res in target_ca:

                start_pos = start_pdb.positions[atom.index].value_in_unit(nanometer)
                target_pos = target_ca[res].value_in_unit(nanometer)
                #print("general", res, start_pos, target_pos, flush=True)
                #print("types", type(start_pos), type(target_pos), flush=True)
                #print("length", len(target_ca), flush=True)
                params = list(start_pos) + list(target_pos)
                #print(params, flush=True)
                dist = np.linalg.norm(start_pos - target_pos)
                print(res, dist, flush=True)
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

    with open(f"{output_prefix}_new_{target_pdb_file}","w") as f:
        PDBFile.writeFile(simulation.topology, state.getPositions(), f)

if __name__ == '__main__':
    #s = time.time()
    args = parser_arg.parse_args()
    #aa_structures = args.aa_structures
    output_prefix = args.output_prefix
    map_small = args.mapping_smaller
    map_big = args.mapping_bigger
    #for input_pdb in inputs:
    #    minimize(input_pdb, integrator, platform, output_prefix)
    #print((time.time()-s))
    map_small_arr = read_mapping(map_small)
    map_big_arr = read_mapping(map_big)
    summary = map_to_map(map_small_arr, map_big_arr)
    summary.to_csv("summary.csv", sep="\t")
    aa_structures = [f"cluster_10/cluster_{i}_isolde_{i}_new_align.pdb" for i in range(10)]
    target_pdbs = [f"cluster_all_new_{i}_align.pdb" for i in range(100)]
    for i, pdb in enumerate(target_pdbs):
        #print(pdb)
        aa_struct_idx = summary.T.idxmax()[i]
        aa_structure = aa_structures[aa_struct_idx]
        #print(aa_structure)
        aa_move_to_ca(aa_structure, pdb, output_prefix) # aa_sturcture,pdb
        print(f"Finished for {pdb}.", flush=True)
