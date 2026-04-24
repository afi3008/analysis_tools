import MDAnalysis as mda
import numpy as np
import pandas as pd
from sklearn.cluster import MiniBatchKMeans
from MDAnalysis.analysis import align
import warnings
import argparse
import re
import os
import glob
from pathlib import Path

parser_arg = argparse.ArgumentParser()
parser_arg.add_argument('--pdb_0', type=str, help='Path to strucutre_z_0.pdb, or alternatively a suitable topology.')
parser_arg.add_argument('--path_to_pdbs', type=str, help='Path to the pdbs to cluster.') # --path_structures
parser_arg.add_argument('--align', action=argparse.BooleanOptionalAction, default=False, help='Align structures first according to RMSD.')
parser_arg.add_argument("--output_path", type=str, required=True, help="path to the pdb file containing the combined pdb")
parser_arg.add_argument('--sele', type=str, default='protein', help="AtomSelection for loading structures. Default is protein.")


warnings.filterwarnings("ignore",
    category=UserWarning,
    module="MDAnalysis.coordinates.TRJ")


def stream_aligned(path_pdb, selection, align=False, pattern="structure_z_*.pdb", chunk_size=1000):
    files = sorted(
        Path(path_pdb).glob(pattern),
        key=lambda p: int(p.stem.split("_")[-1])
    )
    # initialize once
    u = mda.Universe(str(files[0]))
    ref = mda.Universe(str(files[0]))
    atoms = u.select_atoms(selection)
    ref_atoms = ref.select_atoms(selection)

    chunk = []
    
    for p in files:
        u.load_new(str(p))  # load next "frame"
        if align:
            align.alignto(u, ref, select=selection)
        chunk.append(atoms.positions.copy())

        if len(chunk) == chunk_size:
            yield np.array(chunk)
            chunk = []
    if chunk:
        yield np.array(chunk)


if __name__ == "__main__":
    args = parser_arg.parse_args()
    output_path = args.output_path  
    input_top = args.pdb_0 
    struct_path = args.path_to_pdbs
    atomSelection = args.sele
    align = args.align

    kmeans = MiniBatchKMeans(
    n_clusters=4,
    batch_size=1000)

    for chunk in stream_aligned(struct_path, atomSelection, align):
        X = chunk.reshape(chunk.shape[0], -1)
        kmeans.partial_fit(X)

    labels_all = []

    for chunk in stream_aligned(struct_path, atomSelection, align):
        X = chunk.reshape(chunk.shape[0], -1)
        labels = kmeans.predict(X)
        labels_all.append(labels)

    labels_all = np.concatenate(labels_all)

    centroids = kmeans.cluster_centers_

    for i, center in enumerate(centroids):
        coords = center.reshape(-1, 3)

        u = mda.Universe(input_top)
        atoms = u.select_atoms("name CA")
        atoms.positions = coords

        atoms.write(f"cluster_{i}_centroid.pdb")

# Population of clusters
    prob = (
    pd.Series(labels_all)
    .value_counts(normalize=True)
    .sort_index()
    .rename_axis("cluster")
    .reset_index(name="fraction"))
    prob.to_csv("summary_cluster.csv", sep='\t', index=False)

    df = pd.DataFrame({
    "frame": np.arange(len(labels_all)),
    "cluster": labels_all})

    df.to_csv("frame_cluster_map.csv", sep='\t', index=False)
