import numpy as np
import pandas as pd
from typing import List

import os

def rewrite(
    source_dir: str, 
    dest_dir: str,
    idx_ip3_dir_name: int, 
    ip3_dir_names: List[str],
    idx_sample_dir_name: int,
    sample_dir_names: List[str]
    ):

    if idx_sample_dir_name % 100 == 0:
        print("%d / %d - %d / %d" % (
            idx_ip3_dir_name,
            len(ip3_dir_names),
            idx_sample_dir_name,
            len(sample_dir_names)
            ))

    ip3_dir_name = ip3_dir_names[idx_ip3_dir_name]
    sample_dir_name = sample_dir_names[idx_sample_dir_name]
    sample_dir = os.path.join(source_dir,ip3_dir_name,sample_dir_name)
    file_names = os.listdir(sample_dir)

    # Import
    data = {}
    for file_name in file_names:
        base_name = os.path.basename(file_name)
        base_name = os.path.splitext(base_name)[0]

        file_name_full = os.path.join(sample_dir, file_name)
        data[base_name] = pd.read_csv(
            file_name_full,
            index_col=0, 
            delimiter=' ', 
            header=None, 
            names=['t',base_name]
            )

    # Single ordered data frame
    keys = sorted(list(data.keys()))
    data_write = data[keys[0]]
    for key in keys[1:]:
        data_write = data_write.join(data[key])

    # Write
    ip3_dir_name = ip3_dir_names[idx_ip3_dir_name]
    dest_fname = os.path.join(dest_dir,ip3_dir_name,sample_dir_name)+".txt"
    data_write.to_csv(dest_fname, sep=" ")

def rewrite_samples(
    source_dir: str, 
    dest_dir: str, 
    idx_ip3_dir_name: int, 
    ip3_dir_names: List[str]
    ):
    
    # Ensure ip3 dir exists
    ip3_dir_name = ip3_dir_names[idx_ip3_dir_name]
    dest_ip3_dir = os.path.join(dest_dir,ip3_dir_name)
    if not os.path.isdir(dest_ip3_dir):
        os.mkdir(dest_ip3_dir)

    # Dir names
    source_ip3_dir = os.path.join(source_dir,ip3_dir_name)
    sample_dir_names = sorted(os.listdir(source_ip3_dir))

    for idx_sample_dir_name in range(0,len(sample_dir_names)):
        rewrite(
            source_dir=source_dir, 
            dest_dir=dest_dir,
            idx_ip3_dir_name=idx_ip3_dir_name, 
            ip3_dir_names=ip3_dir_names,
            idx_sample_dir_name=idx_sample_dir_name,
            sample_dir_names=sample_dir_names
            )

def rewrite_ip3rs(source_dir: str, dest_dir: str):

    # Make destination dir if it doesnt exist
    if not os.path.isdir(dest_dir):
        os.mkdir(dest_dir)

    ip3_dir_names = sorted(os.listdir(source_dir))

    for idx_ip3_dir_names in range(0,len(ip3_dir_names)):
        rewrite_samples(
            source_dir=source_dir, 
            dest_dir=dest_dir, 
            idx_ip3_dir_name=idx_ip3_dir_names, 
            ip3_dir_names=ip3_dir_names
            )

if __name__ == "__main__":

    source_dir = "../data_tau_leaping/vol_exp_14/ip3r_01000/"
    dest_dir = "../data_tau_leaping_rw/vol_exp_14/ip3r_01000/"

    rewrite_ip3rs(
        source_dir=source_dir,
        dest_dir=dest_dir
    )