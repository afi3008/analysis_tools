import numpy as np
import pandas as pd
import starfile
import argparse


parser_arg = argparse.ArgumentParser()
parser_arg.add_argument('--path_starfile', type=str, required=True)
parser_arg.add_argument('--mapping_path', type=str)
parser_arg.add_argument('--prefix_output_starfile', type=str)
parser_arg.add_argument('--stride', type=int, default=1)
args = parser_arg.parse_args()
path_starfile = args.path_starfile
mapping_path = args.mapping_path
prefix_output_starfile = args.prefix_output_starfile
stride = args.stride

def read_mapping(mapping_path):
    """
    Read mapping from clustering (which particle belongs to which cluster).
    """
    cluster_mapping = np.genfromtxt(mapping_path, delimiter=None, dtype=int)
    mapping = cluster_mapping[:, 1]
    return mapping

def read_starfile(path_starfile):
    """
    Read starfile.
    """
    df = starfile.read(path_starfile)
    return df

def stride_starfile(df_starfile, stride):
    """
    When dataset to large to cluster as whole, clustering only performed on every 2nd or 3rd structure, called striding. Adapt starfile according to that before filtering.
    :param stride: Every n-th structure that was used for clustering.
    """
    df_stride = df_starfile.copy()
    particles = df_stride['particles']
    particles = particles[particles.index % stride == 0]
    df_stride['particles'] = particles
    return df_stride

def cluster_starfile(df_starfile, cluster):
    """
    Filter starfile according to cluster mapping. Remove all particles that are not in the sought after cluster.
    :param cluster: The sought after cluster the dataset gets filtered upon.  
    """
    df_cluster = df_starfile.copy()
    particles = df_cluster['particles']  # Access particles df of starfile.
    particles = particles[particles['cluster'] == cluster]  # Filtering.
    particles = particles.drop(columns='cluster')  # Drop cluster column.
    df_cluster['particles'] = particles
    return df_cluster

if __name__ == '__main__':
    mapping = read_mapping(mapping_path)
    df_star = read_starfile(path_starfile)
    df_star = stride_starfile(df_star, stride)
    df_star['particles']['cluster']= mapping  # Add mapping as column to the pd.df.
    
    for cluster in range(df_star['particles']['cluster'].max()+1):  # Only works if numbering of clusters is sequential.
        df_cluster = cluster_starfile(df_star, cluster)
        starfile.write(df_cluster, f"{prefix_output_starfile}_cluster_{cluster}.star")
