import argparse
import os

import numpy as np
import pandas as pd
from modify_cluster_output_file import *


def get_lc_reads(output_lst : list):
    """
    Modify sample A and sample B 's clubcpg output for extracting lc_reads
    :param output_lst:  a list of path names
    :return: a modified dataframe?
    """

    club_combined = pd.DataFrame()
    file_idx = 1
    for output_path in output_lst:
        clubcpg = pd.read_csv(output_path)

        clubcpg = cluster_output_change_bin_name(cluster_df=clubcpg)
        clubcpg = cluster_output_separate_AB(cluster_df=clubcpg)
        clubcpg = cluster_output_add_start_end(cluster_df=clubcpg)
        clubcpg["origin"] = file_idx # each file index represent a sample

        club_combined = pd.concat([club_combined, clubcpg])
        print("done")
        file_idx += 1

    summaryA = club_combined.groupby(["bin_id", "origin"], as_index=False)["A"].sum()
    summaryB = club_combined.groupby(["bin_id", "origin"], as_index=False)["B"].sum()
    summaryA = summaryA.groupby("bin_id", as_index=False)["A"].min()
    summaryB = summaryB.groupby("bin_id", as_index=False)["B"].min()
    summaryAB = pd.concat([summaryA, summaryB["B"]], axis=1)
    summaryAB["lc_min"] = summaryAB[["A", "B"]].min(axis=1)
    summaryAB = summaryAB[["bin_id", "lc_min"]]
    return summaryAB

### TEST ###
path1 = "/Users/david/Sphere_files/Downsample replicate/Clubcpg_re_run_July2021/sex_p35_neuron/male_p35_neuron.bam.chr19_cluster_results.csv"
path2 = "/Users/david/Sphere_files/Downsample replicate/Clubcpg_re_run_July2021/sex_p12_neuron/male_p12_neuron.bam.chr19_cluster_results.csv"


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # parser.add_argument("-A", help="absolute path to the first clubcpg cluster output", default=None)
    # parser.add_argument("-B", help="absolute path to the second clubcpg cluster output", default=None)
    parser.add_argument("--samples", help="absolute path to the  clubcpg cluster outputs", default=None, nargs="+")
    parser.add_argument("-o", "--output", help="folder to save imputed coverage data", default=None)
    args = parser.parse_args()

    if not args.output:
        output_folder = os.path.dirname(args.A + "/../DEC7")
    else:
        output_folder = args.output
    try:
        os.mkdir(output_folder)
    except FileExistsError:
        print("Output folder already exists... no need to create it...")

    lc_reads = get_lc_reads(args.samples)
    lc_reads.to_csv(args.output + "/lc_reads1.csv")
    print(lc_reads)
    print("done")


