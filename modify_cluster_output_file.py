import pandas as pd
import numpy as np
import argparse
import os
from collections import defaultdict
def cluster_output_change_bin_name(cluster_df:pd.DataFrame):
    cluster_df.rename(columns={'bin': "bin_id"}, inplace=True)
    return cluster_df

def cluster_output_separate_AB(cluster_df: pd.DataFrame):
    """
    get a cluster_output.csv 's dataframe, do the following changes:
    * split input lable from AB to 1 col of A and 1 col of B
    :param cluster_df: a cluster_output.csv 's dataframe
    :return: the modified dataframe??
    """
    # get all rows with label AB
    temp_AB = cluster_df.loc[cluster_df["input_label"] == "AB"].copy()
    # separate class split into two rows
    temptemp = temp_AB["class_split"].str.split(";", expand=True)
    temp_AB["A"] = temptemp[0]
    temp_AB.loc[:, "B"] = temptemp[1]
    # get all rows with label A
    temp_A = cluster_df.loc[cluster_df["input_label"] == "A"].copy()
    temp_A.loc[:, "A"] = temp_A.loc[:, "class_split"]
    temp_A.loc[:, "B"] = "B=0"
    # get all rows with label B
    temp_B = cluster_df.loc[cluster_df["input_label"] == "B"].copy()
    temp_B.loc[:, "A"] = "A=0"
    temp_B.loc[:, "B"] = temp_B.loc[:, "class_split"]
    # combine all rows of label A,B,AB
    temp_all_labels = pd.concat([temp_AB, temp_A, temp_B])
    # get rid of the "A=", and "B=" in the two cols
    temp_all_labels["A"] = temp_all_labels["A"].str.replace("A=", "", regex=True)
    temp_all_labels["B"] = temp_all_labels["B"].str.replace("B=", "", regex=True)
    # turn vals in col A, B into numbers;
    temp_all_labels["A"] = temp_all_labels["A"].astype(int)
    temp_all_labels["B"] = temp_all_labels["B"].astype(int)

    return temp_all_labels.reset_index(drop=True)


def cluster_output_add_start_end(cluster_df: pd.DataFrame):
    """
    separate the 'bin' column into three cols: chr, start, end
    :param cluster_df: the clubcpg output file ( other columns might be modified, such as in the case if separate_AB is
    called
    :return: the modified dataframe
    """
    temp = cluster_df["bin_id"].str.split("_", expand=True).copy()
    cluster_df["chr"] = temp[0]
    cluster_df["end"] = temp[1].astype(int)
    cluster_df["start"] = cluster_df["end"] - 99
    # NOTE: Here -99 because in harry's paper it was 100 bin size.
    # Fixme: how to allow modifying this ↑?
    # reorder columns:
    cluster_df = cluster_df[['bin_id', 'chr', 'start', 'end', 'input_label', 'methylation', 'class_label', 'read_number',
                             'cpg_number', 'cpg_pattern', 'class_split', 'A', 'B']]
    return cluster_df


def cluster_output_replace_label(cluster_df: pd.DataFrame, patterns: pd.DataFrame):
    """
    replace the original clubcpg class labels with consistent class labels across different bins. In Harry's script:
    "# Anthony's class_label variable isn't so great (it seems to determine them either per-bin or per-analysis, but I'd
     rather have global ones that apply to all my datasets)"
    :param cluster_df: the clubcpg output file ( other columns might be modified, such as in the case if separate_AB is
    called
           patterns: the new patterns dataframe from lc_reads files, consistent among bins
    :return: the modified dataframe
    """
    patterns = patterns[["cpg_pattern", "class_label"]]
    cluster_df = cluster_df.loc[:, cluster_df.columns != "class_label"].merge(patterns,
                                                                              how="left", on="cpg_pattern")
    # reorder columns:
    cluster_df = cluster_df[['bin_id', 'chr', 'start', 'end', 'input_label', 'cpg_number', 'class_label', 'methylation',
                             'read_number', 'cpg_pattern', 'class_split', 'A', 'B']]
    return cluster_df

def cluster_output_add_lcreads(cluster_df: pd.DataFrame, lc_reads:pd.DataFrame):
    """
    add the lc_read of each bin to the dataframe
    :param cluster_df: the clubcpg output file ( other columns might be modified, such as in the case if separate_AB is
    called
    :param lc_reads: the lc_read dataframe at each bin
    :return: the modified dataframe
    """
    lc_reads = lc_reads[["bin_id","lc_min"]]
    lc_reads = lc_reads.rename(columns={"lc_min": "lc_sum"})

    cluster_df = cluster_df.merge(lc_reads,how="left", on="bin_id")
    return cluster_df

def cluster_output_add_v1(cluster_df:pd.DataFrame):
    """
    sort the df by end
    add V1 variable to track process: each V1 value is a percentage, each bin has one V1 value, regardless of epiallele;
    :param cluster_df: the clubcpg output file ( other columns might be modified, such as in the case if separate_AB is
    called
    :return: the modified dataframe
    """
    # sort cluster_df based on end position
    cluster_df.sort_values("bin_id", inplace=True)
    cluster_df.reset_index(drop=True,inplace=True)
    # group bins and add a V1 val for each bin
    bins = cluster_df["bin_id"].unique()
    bins = pd.DataFrame(bins)
    bins.rename(columns={0:"bin_id"},inplace=True)
    bins["V1"] = bins.index+1
    cluster_df = cluster_df.merge(bins,how="left", on="bin_id")
    cluster_df["V1"] = cluster_df["V1"].astype(float)/bins.shape[0] * 100
    # print(bins)
    return cluster_df

# ### TEST ###
# TEST_CSV1 = "/Users/david/Sphere_files/Downsample/Clubcpg_re_run_July2021/sex_p35_neuron/male_p35_neuron.bam.chr19_cluster_results.csv"
# TEST_CSV2 = "/Users/david/Sphere_files/Downsample/Clubcpg_re_run_July2021/sex_p12_neuron/male_p12_neuron.bam.chr19_cluster_results.csv"
# TEST_CSV_sample = "/Users/david/Sphere_files/Downsample replicate/CluBCpG demos/raw_data/sample clubcpg output.csv"
# sample_df = pd.read_csv(TEST_CSV_sample)
# modified_df = cluster_output_separate_AB(sample_df)
# # print(modified_df)
# modified_df2 = cluster_output_add_start_end(modified_df)
# # print(modified_df2)
# patterns = pd.read_csv("/Users/david/Sphere_files/Downsample replicate/CluBCpG demos/output_csv/cluster_patterns.csv")
# modified_df3 = cluster_output_replace_label(modified_df2,patterns)
# # print(modified_df3)
# lc_reads = pd.read_csv("/Users/david/Sphere_files/Downsample replicate/CluBCpG demos/output_csv/lowest common read depths.csv")
# modified_df4 = cluster_output_add_lcreads(modified_df3,lc_reads)
# # print(modified_df4)
# modified_df5 = cluster_output_add_v1(modified_df4)
# # print(modified_df5)
# ### TEST END ###
