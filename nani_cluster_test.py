import MDAnalysis as mda
import numpy as np
import pandas as pd
from sklearn.cluster import MiniBatchKMeans
from MDAnalysis.analysis import align
import warnings
import argparse


parser_arg = argparse.ArgumentParser()
parser_arg.add_argument('--pdb_0', type=str, help='Path to strucutre_z_0.pdb, or alternatively a suitable topology.')
parser_arg.add_argument('--path_to_pdbs', type=str, help='Path to the pdbs to cluster.')
parser_arg.add_argument('--align', action=argparse.BooleanOptionalAction, default=False, help='Align structures first according to RMSD.')

warnings.filterwarnings("ignore",
    category=UserWarning,
    module="MDAnalysis.coordinates.TRJ")


def stream_traj_numpy(top, traj, selection, chunk_size=1000):
    u = mda.Universe(top, traj)
    atoms = u.select_atoms(selection)

    chunk = []
    for ts in u.trajectory:
        coords = atoms.positions.copy()  # (n_atoms, 3)
        chunk.append(coords)

        if len(chunk) == chunk_size:
            yield np.array(chunk)  # (chunk_size, n_atoms, 3)
            chunk = []

    if chunk:
        yield np.array(chunk)


def stream_aligned(top, traj, selection, chunk_size=1000):
    u = mda.Universe(top, traj)
    ref = mda.Universe(top, traj)  # first frame as reference
    align.AlignTraj(u, ref, select=selection, in_memory=False).run()
    atoms = u.select_atoms(selection)
    chunk = []
    for ts in u.trajectory:
        chunk.append(atoms.positions.copy())
        if len(chunk) == chunk_size:
            yield np.array(chunk)
            chunk = []
    if chunk:
        yield np.array(chunk)


if __name__ == '__main__':
    args = parser_arg.parse_args()
    input_top = arg.pdb_0 
    input_traj = "full.nc" 
    atomSelection = "protein"

kmeans = MiniBatchKMeans(
    n_clusters=4,
    batch_size=1000)

for chunk in stream_traj_numpy(input_top, input_traj, atomSelection):
    X = chunk.reshape(chunk.shape[0], -1)
    kmeans.partial_fit(X)

labels_all = []

for chunk in stream_traj_numpy(input_top, input_traj, atomSelection):
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
    .reset_index(name="fraction")
)
prob.to_csv("summary_cluster.csv", sep='\t', index=False)

df = pd.DataFrame({
    "frame": np.arange(len(labels_all)),
    "cluster": labels_all
})

df.to_csv("frame_cluster_map.csv", sep='\t', index=False)
