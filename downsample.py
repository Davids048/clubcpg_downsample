import itertools
import math
from functools import partial

# import dask.dataframe as dd
import pandas as pd
import numpy as np
import time
import scipy as sp
from modify_cluster_output_file import *
from multiprocessing import Pool, cpu_count
# from dask.distributed import LocalCluster, Client
import multiprocessing
from typing import Callable, Tuple, Union
from datetime import date


start_time = time.time()
print("start",start_time)
progress = 1
def fisher_on_df(combined_bins: pd.DataFrame, read_target: int):
    """
    given the combined bins (pivot_longer applied), return the -log p val from each row, representing if there
    is a significant difference between library A and library B
    :param combined_bins: combined bins (pivot_longer applied),
    :param read_target: lc_read depth at that bin. / amount of labels that got chosen in each sampling iteration.
    :return: a df, with p values of each row.
    """
    # TODO: implement this
    fisher_df = pd.DataFrame(columns=["A", "B"])
    fisher_df["A"] = combined_bins["A"]
    # fisher_df["read_target"] = read_target
    fisher_df["B"] = combined_bins["B"]
    # CHECK: if this is the correct fisher val, and if it is the same as r script
    # row = fisher_df.iloc[9]
    # print("this array\n", np.array([[row["A"], row["B"]], row["read_target"]], dtype=object))
    # print("pval\n",sp.stats.fisher_exact(np.array([[7,12],[14,5]], dtype=object),alternative='two-sided')[1])
    # CheckEnd
    # Fixed:  error: fisher_df["p_val"] = fisher_df.apply(lambda row: sp.stats.fisher_exact(
    #  raise ValueError("All values in `table` must be nonnegative.")
    # ValueError: All values in `table` must be nonnegative.
    ######debug#####
    # print(fisher_df)
    # for row_idx in  range(fisher_df.shape[0]):
    #     x = fisher_df.loc[row_idx]
    #     print(x)
    #     print("target",read_target)
    #     print(np.array([[x["A"], read_target - x["A"]], [x["B"], read_target - x["B"]]]))

    #####debugend#####
    fisher_df["p_val"] = fisher_df.apply(lambda row: sp.stats.fisher_exact(
        np.array([[row["A"], read_target - row["A"]], [row["B"], read_target - row["B"]]], dtype=object),
        alternative='two-sided')[1], axis=1)
    fisher_df["p_val"] = fisher_df.apply(lambda row: (-1) * math.log10(row["p_val"]), axis=1)
    # print(fisher_df)
    # TODO: fisher is not right, modify
    return fisher_df


def bin_resample(bin_name, df: pd.DataFrame):
    # original header: def bin_resample(df: pd.DataFrame, n: int):
    n = 100
    """
    resample for each bin
    :param df: a dataframe of methylation information at a bin. if there are n epialleles, there should be n rows.
    :param n: the number of resampling to do
    :return: a new dataframe with the same information as the final csv for down sampling.
    """
    # if cluster has one epi-allele, can skip all things
    # print("progress?", df["V1"].copy().reset_index(drop=True)[0])
    global progress
    # print(progress)
    progress +=1
    # print(df["bin_id"])
    # print(df["lc_sum"])
    # print(df["A"])
    # print(df["B"])
    # print(PROGRESS/Size*100)
    # PROGRESS +=1
    if df.shape[0] == 1:
        # CHECK: if the output is same as row != 1
        lc_read = df["lc_sum"].copy().reset_index(drop=True)[0]
        data = {"class_label": [max(df["class_label"])],
                "A_count_mean": [lc_read],
                "B_count_mean": [lc_read],
                "A_norm_mean": [100],
                "B_norm_mean": [100],
                "delta_mean": [0],
                "A_count_sd": [0],
                "B_count_sd": [0],
                "A_norm_sd": [0],
                "B_norm_sd": [0],
                "p_val_mean": [1],
                "sig_pct": [0]
                }
        res =  pd.DataFrame(data)
        res["bin_id"] = bin_name
        return  res
    else:
        # if cluster has multiple epi-alleles:
        # initialize
        # read_target = df["lc_sum"][1]
        # print(df)
        # print(type(df["lc_sum"]))
        read_target = df["lc_sum"].reset_index(drop=True)[0]
        class_max = max(df["class_label"])
        # create the array of elements to sample from
        reads_A = np.array(df["class_label"])
        reads_A = np.repeat(reads_A, df["A"], axis=0)
        reads_B = np.array(df["class_label"])
        reads_B = np.repeat(reads_B, df["B"], axis=0)
        labels = np.sort(df["class_label"].unique())

        # print(reads_A)

        def sample_from_reads(reads: np.ndarray, target: int, size: int):
            """
            Given an array of read lables, sample targeted_number of reads from the array n times
            put each sample into a dataframe of #_class_max cols, and n rows, representing the n samples.
            :return: a dataframe of sampled labels' counts
            """
            lst_sampls = [np.random.choice(reads, size=target, replace=False) for i in range(size)]
            # print(lst_sampls)

            # for each sample, put into a row of the output dataframe. To ensure consistency, the dataframe
            # will have all possible lables in this bin, by extracting all unique labels from df["class_labels"]
            labels = np.sort(df["class_label"].unique())
            table = pd.DataFrame(index=labels)
            # TODO: collect distribution of labels in each row of lst_samples, put into table.
            # for each elem (sampled labels) in lst_samples:
            for sample in lst_sampls:
                counts = pd.Categorical(sample, categories=labels, ordered=True).value_counts()
                table = pd.concat([table, counts], axis=1)
            # # cat = pd.Categorical(lst_sampls[1], categories=[i for i in range(class_max+1)], ordered=True)
            # cat = pd.Categorical(lst_sampls[1], categories=labels, ordered=True)
            # counts = cat.value_counts()
            # print("counts\n", counts)
            # table = pd.concat([table,counts], axis = 1)
            #
            # cat = pd.Categorical(lst_sampls[4], categories=labels, ordered=True)
            # counts = cat.value_counts()
            # print("counts\n", counts)
            # table = pd.concat([table, counts], axis=1)
            # print(table.T.reset_index(drop=True))
            table = table.T.reset_index(drop=True)
            table.index += 1
            table["iteration"] = table.index
            return table

        def helper_create_corresponding_table(reads: np.ndarray, size: int):
            """
            create a df that has the same number of rows as the sampling df. But without doing any sampling,
            cuz it's just a waste a time to sample #target things from #target things.
            :param reads: the reads to sample from
            :param size: the number of samples to do, should be same as the size parameter as sample_from_reads
            :return: a dataframe.
            """
            labels = np.sort(df["class_label"].unique())
            row = pd.Categorical(reads, categories=labels, ordered=True).value_counts()
            table = pd.DataFrame(row, index=labels)
            # print(table.T)
            table = pd.DataFrame(np.repeat(table.values, size, axis=1))
            table = table.T.reset_index(drop=True)
            table.columns = labels
            # print("tableows\n", table)
            table.index += 1
            table["iteration"] = table.index
            return table

        # make the two tables
        if sum(df["A"]) < read_target or sum(df["B"]) < read_target:
            if sum(df["A"]) < sum(df["B"]):
                # print(1)
                table_A = helper_create_corresponding_table(reads_A, n)
                # print(table_A)
                table_B = sample_from_reads(reads_B, read_target, n)
                # print(table_B)
            elif sum(df["A"]) > sum(df["B"]):
                # print(2)
                table_B = helper_create_corresponding_table(reads_B, n)
                # print(table_A)
                table_A = sample_from_reads(reads_A, read_target, n)
                # print(table_B)
            else:
                # print(3)
                table_A = helper_create_corresponding_table(reads_A, n)
                table_B = helper_create_corresponding_table(reads_B, n)
        else:
            # print(4)
            table_A = sample_from_reads(reads_A, read_target, n)
            table_B = sample_from_reads(reads_B, read_target, n)

        # turn tables from wide to long
        bins_A = pd.melt(table_A, id_vars=["iteration"], value_vars=labels, var_name=["class_label"]).sort_values( \
            by=["iteration", "class_label"]).reset_index(drop=True).rename(columns={"value": "A"})
        # print("bins_A\n", bins_A)
        bins_B = pd.melt(table_B, id_vars=["iteration"], value_vars=labels, var_name=["class_label"]).sort_values( \
            by=["iteration", "class_label"]).reset_index(drop=True).rename(columns={"value": "B"})
        # print("bins_B\n", bins_B)
        combined_bins = pd.concat([bins_A, bins_B["B"]], axis=1)
        # print(combined_bins)

        # redo clustering for each sampling
        combined_bins["cluster_label"] = "AB"
        combined_bins.loc[(combined_bins["A"] >= 4) & (combined_bins["B"] == 0), ["cluster_label"]] = "A"
        combined_bins.loc[(combined_bins["B"] >= 4) & (combined_bins["A"] == 0), ["cluster_label"]] = "B"
        # print('hi\n',combined_bins)

        # print(combined_bins)
        # print("read target", read_target)
        # calculate enrichment
        fisher_df = fisher_on_df(combined_bins, read_target)
        # print(fisher_df)
        combined_bins = pd.concat([combined_bins, fisher_df["p_val"]], axis=1)
        # print(combined_bins)
        # give a discrete enrichment identifier
        combined_bins["is_sig"] = 0
        # combined_bins["is_sig"][combined_bins["p_val"] > (-1)*math.log10(0.05)] = 1
        combined_bins.loc[combined_bins["p_val"] > (-1) * math.log10(0.05), 'is_sig'] = 1
        # print(combined_bins)

        # summarize stats
        combined_sum = pd.DataFrame()
        combined_sum["class_label"] = labels
        combined_sum.index = labels
        grouped_df = combined_bins.groupby(["class_label"], as_index=True)
        # CHECK if calculation is right
        # combined_sum["A_count_mean"] = combined_bins.groupby(["class_label"], as_index=True).agg({"A":"mean"})["A"]
        combined_sum["A_count_mean"] = grouped_df["A"].mean()
        combined_sum["B_count_mean"] = grouped_df["B"].mean()
        combined_sum["A_norm_mean"] = grouped_df["A"].mean() / read_target
        combined_sum["B_norm_mean"] = grouped_df["B"].mean() / read_target * 100
        combined_sum["delta_mean"] = (grouped_df["A"].mean() - grouped_df["B"].mean()) / read_target * 100
        combined_sum["A_count_sd"] = grouped_df["A"].std()
        combined_sum["B_count_sd"] = grouped_df["B"].std()
        combined_sum["A_norm_sd"] = grouped_df["A"].std() / read_target * 100
        combined_sum["B_norm_sd"] = grouped_df["B"].std() / read_target * 100
        # print("watereve\n" , ((grouped_df["A"].std() + grouped_df["B"].std()) / read_target) * 100)
        combined_sum["delta_sd"] = ((grouped_df["A"].std() + grouped_df["B"].std()) / read_target) * 100
        combined_sum["delta_sd"] = ((grouped_df["A"].std() + grouped_df["B"].std()) / read_target) * 100
        combined_sum["p_val_mean"] = grouped_df["p_val"].mean()
        # print(grouped_df["is_sig"].describe())
        combined_sum["sig_pct"] = grouped_df["is_sig"].mean() * 100

        #
        #
        # print("summary norm mean\n", combined_sum["delta_mean"])
        # print("summary\n", combined_sum)
        # print(combined_bins.groupby(["class_label"], as_index=True)["B"].describe())
        # print("combinedsum\n", combined_bins)
        combined_sum["bin_id"] = bin_name
        return combined_sum


def downsampel_prep(patterns_path: str, lc_reads_path: str, clubcpg_path_A: str, clubcpg_path_B: str, sampleA:str, sampleB:str):
    """
    Process the input class-patterns, lowest-common-read depth, clubcpg output for bin_resample function;
    :param patterns_path: path for the class label file
    :param lc_reads_path: path for the file storing lowest common read depth of each bin
    :param clubcpg_path: path for file that stores clubcpg clustering output
    :param sampleA: which sample in cluster output A to choose from, allowed inputs are either A or B
    :param sampleB: which sample in cluster output B to choose from, allowed inputs are either A or B
    :return: the final result after resampling
    """
    patterns = pd.read_csv(patterns_path)
    lc_reads = pd.read_csv(lc_reads_path)
    # print("lc shape",lc_reads.shape)
    # print(lc_reads)
    club_A = pd.read_csv(clubcpg_path_A)
    club_B = pd.read_csv(clubcpg_path_B)

    # print("A shape", club_A.shape)
    # print("B shape", club_B.shape)

    club_combined = pd.DataFrame()

    # align bin naming: done in method already
    # clubcpg.rename(columns={"bin": "bin_id"}, inplace=True)

    file_idx = 1
    for club in [club_A, club_B]:
        # select bins that meets lc_reads bins;
        # CHECK: harry's script chose only bins within the lc_reads thing. Is that really necessary?
        club = cluster_output_change_bin_name(cluster_df=club)
        club = club[club["bin_id"].isin(lc_reads["bin_id"])]
        club = cluster_output_separate_AB(club)
        club = cluster_output_add_start_end(cluster_df=club)
        club = cluster_output_replace_label(cluster_df=club, patterns=patterns)
        club["origin"] = file_idx
        # this column denote where the file is from， 1 means the file from -A input file, 2 means
        # from -B file. In the current analysis, we are looking for male, p35 vs p12, so -A is p35, -B is p12
        # (or vice versa)
        # in each dataframe, male is sample A, so we take the A colomn, (if want female, use B coloumn)
        file_idx += 1
        # print("selected", club.shape)
        club_combined = pd.concat([club_combined, club])

    club_index = create_idx_file(club_combined)
    ## Fixme: assume we are doing male, (will need to generalize later),
    #   i.e. doing male, so select all bins with
    club_a = club_combined[club_combined["origin"] == 1]
    # club_a = club_a[["bin_id", "class_label", "A"]]
    # in each dataframe, male is sample A, so we take the A column, (if want female, use B column)
    club_a = club_a[["bin_id", "class_label", sampleA]]
    club_a = club_a.rename(columns={sampleA: "B"})
    club_a = club_a.drop_duplicates()
    # club_a.rename(columns={"origin":"A"})
    club_b = club_combined[club_combined["origin"] == 2]
    club_b = club_b[["bin_id", "class_label", sampleB]]
    club_b = club_b.rename(columns={sampleB: "B"})  # in the new df, it will be the second thing to compare
    club_b = club_b.drop_duplicates()
    # print(club_b)
    # club_b.rename(columns={"origin": "B"})
    # print(club_b)
    club_merge = club_a.merge(club_b, how="outer", on=["bin_id", "class_label"]).fillna(0)
    club_merge = cluster_output_add_lcreads(cluster_df=club_merge, lc_reads=lc_reads)

    club_merge = cluster_output_add_v1(cluster_df=club_merge)
    print("done prep func")
    return club_index, club_merge


def downsampe_prep_single_file(patterns_path: str, lc_reads_path: str, clubcpg_path: str): # -> clubidx, clubcpg
    """
    Process the input class-patterns, lowest-common-read depth, clubcpg output for bin_resample function;
    :param patterns_path: path for the class label file
    :param lc_reads_path: path for the file storing lowest common read depth of each bin
    :param clubcpg_path: path for file that stores clubcpg clustering output
    :return: the final result after resampling
    """
    patterns = pd.read_csv(patterns_path)
    lc_reads = pd.read_csv(lc_reads_path)
    # print(lc_reads)
    clubcpg = pd.read_csv(clubcpg_path)

    # align bin naming: done in method already
    # clubcpg.rename(columns={"bin": "bin_id"}, inplace=True)

    # select bins that meets lc_reads bins;
    # clubcpg = clubcpg[clubcpg["bin_id"] in lc_reads["bin_id"]]
    # print(lc_reads["bin_id"])
    # clubcpg = clubcpg.loc[~clubcpg.bin_id.isin(lc_reads["bin_id"]).dropna(),:]
    # CHECK: harry's script chose only bins within the lc_reads thing. Is that really necessary?
    clubcpg = cluster_output_change_bin_name(clubcpg)
    clubcpg = cluster_output_separate_AB(clubcpg)
    clubcpg = cluster_output_add_start_end(cluster_df=clubcpg)
    clubcpg = cluster_output_replace_label(cluster_df=clubcpg, patterns=patterns)
    clubcpg = cluster_output_add_lcreads(cluster_df=clubcpg, lc_reads=lc_reads)
    clubcpg = cluster_output_add_v1(cluster_df=clubcpg)
    club_idx = create_idx_file(clubcpg)
    return club_idx, clubcpg


def create_idx_file(clubcpg: pd.DataFrame):
    """
    an index file to combine all bins after applying downsample at each bin
    :param clubcpg: the dataframe after doing downsample prep
    :return: a modified dataframe
    """
    club_idx = clubcpg[["bin_id", "chr", "start", "end", "cpg_number", "class_label", "methylation", "cpg_pattern"]]
    club_idx = club_idx.drop_duplicates()
    return club_idx

def sample_each_bin(small_group1):  # TODO: or groupby?
    res_df1 = small_group1.apply(lambda x: bin_resample(x, 100))
    return res_df1

def run_downsample3(clubcpg_df: pd.DataFrame,club_idx:pd.DataFrame):
    useful_part = clubcpg_df[["bin_id", "class_label", "A", "B", "lc_sum", "V1"]]
    bin_groups = useful_part.groupby(by="bin_id", sort=False)
    print("arrived")
    with Pool(cpu_count()) as p:
        args = [(group, 100) for group in bin_groups]
        # ret_list = p.map(bin_resample,[group for name, group in bin_groups], [100]*bin_groups.ngroups)
        ret_list = p.starmap_async(bin_resample, args)
    res = ret_list
    print("finish")
    print(type(res))
    print(res)
    return res

def groupby_parallel(groupby_df: pd.core.groupby.DataFrameGroupBy,
                     func,
                     num_cpus: int,
                     logger: Callable[[str], None]=print) -> pd.DataFrame:
    """Performs a Pandas groupby operation in parallel.
    Example usage:
        import pandas as pd
        df = pd.DataFrame({'A': [0, 1], 'B': [100, 200]})
        df.groupby(df.groupby('A'), lambda row: row['B'].sum())
    Authors: Tamas Nagy and Douglas Myers-Turnbull
    """
    start = time.time()
    logger("\nUsing {} CPUs in parallel...".format(num_cpus))
    with multiprocessing.Pool(num_cpus) as pool:
        queue = multiprocessing.Manager().Queue()
        # result = pool.starmap_async(func, [(name, group) for name, group in groupby_df])
        result = pool.starmap_async(func, [(name,group,) for name, group in groupby_df])
        cycler = itertools.cycle('\|/―')
        while not result.ready():
            # print(queue.qsize()/len(groupby_df))
            logger("Percent complete: {:.0%} {}".format(queue.qsize()/len(groupby_df), next(cycler)), end="\r")
            time.sleep(0.01)
        got = result.get()

    logger("\nProcessed {} rows in {:.1f}s".format(len(got), time.time() - start))
    return pd.concat(got)

def run_downsample(clubcpg_df: pd.DataFrame, club_idx: pd.DataFrame):
    useful_part = clubcpg_df[["bin_id", "class_label", "A", "B", "lc_sum", "V1"]]
    useful_part["bin_id_backup"] = useful_part["bin_id"]
    bin_groups = useful_part.groupby(by="bin_id_backup", sort=False)
    res = groupby_parallel(bin_groups, bin_resample, num_cpus=int(args.ncore))
    print("shape",res)
    res_df = club_idx.merge(res, how="right", on=["bin_id", "class_label"])
    return res_df

def run_downsample4(clubcpg_df: pd.DataFrame, club_idx: pd.DataFrame):
    clust = LocalCluster()
    clt = Client(clust, set_as_default=True)
    useful_part = clubcpg_df[["bin_id", "class_label", "A", "B", "lc_sum", "V1"]]
    dask_useful_part = dd.from_pandas(useful_part, npartitions=2)
    print(dask_useful_part)
    output_df = (dask_useful_part.groupby("bin_id").apply(lambda x:bin_resample(x, 100))).compute(num_workers=1000)

    return output_df

# right, slow, original
def run_downsample2(clubcpg_df: pd.DataFrame, club_idx: pd.DataFrame):
    """
    after preparing input df using downsample_prep, run the downsample thingy on each bin
    :param clubcpg_df: the modified dataframe
    :param club_idx: the idx file used to bind all outputs
    :return: the output dataframe
    """
    useful_part = clubcpg_df[["bin_id", "class_label", "A", "B", "lc_sum", "V1"]]
    # test 2 groups
    # two_opts = ["chr1_185043700"]
    # useful_part = useful_part[useful_part["bin_id"].isin(two_opts)]
    # useful_idx = club_idx[["bin_id","chr","start","end","cpg_number","class_label","methylation","cpg_pattern"]]
    # useful_idx = useful_idx[useful_idx["bin_id"].isin(two_opts)]
    # print(useful_part)
    # DONE: after testing, get rid of the two opts
    bin_groups = useful_part.groupby(by="bin_id", sort=False)
    output_df1 = bin_groups.apply(lambda x: bin_resample(x))
    print(output_df1)
    # print(output_df)
    res_df = club_idx.merge(output_df1, how="right", on=["bin_id", "class_label"])
    # res_df = useful_idx.merge(output_df, how="left", on=["bin_id", "class_label"])
    return res_df


### Tests ###
# bin_multi_allele_data = {"bin_id":["chr1_119077300", "chr1_119077300", "chr1_119077300", "chr1_119077300", "chr1_119077300", "chr1_119077300", "chr1_119077300"],
#             "class_label":[14,13,2,0,11,12,10],
#             "A":[7,2,2,1,3,4,0],
#             "B":[82,4,2,5,4,0,5],
#             "lc_sum":[19,19,19,19,19,19,19],
#             "V1":[4,4,4,4,4,4,4]}
# bin_multi_allele = pd.DataFrame(bin_multi_allele_data)
# print("cool\n",bin_multi_allele)
# print("whoooo\n", bin_resample(bin_multi_allele, 100))
#
# bin_single_allele_data = {"bin_id":["chr1_119077300"],
#                           "class_label":[14],
#                           "A":[7],
#                           "B":[82],
#                           "lc_sum":[19],
#                           "V1":[4]}
# bin_single_allele = pd.DataFrame(bin_single_allele_data)
# # print(bin_resample(bin_single_allele, 2))
### TEST END ###
### TEST 2 ###
## TEST 3 sample file
# if __name__ =="__main__":
#     pattern_path = "/Users/david/Sphere_files/Downsample replicate/CluBCpG demos/output_csv/cluster_patterns.csv"
#     lc_path = "/Users/david/Sphere_files/Downsample replicate/CluBCpG demos/output_csv/lowest common read depths.csv"
#     clubcpg_path = "/Users/david/Sphere_files/Downsample replicate/CluBCpG demos/raw_data/sample clubcpg output.csv"
#     clubcpg = downsampe_prep_single_file(pattern_path, lc_path, clubcpg_path)
#     Size = clubcpg.shape[0]
#     print("finished prep")
#     club_idx = create_idx_file(clubcpg)
#     print(club_idx.columns)
#     print("finished idx")
#     output_df = run_downsample(clubcpg,club_idx)
#     print(output_df)
#     print("time:", time.time()- start_time, "s")
### TEST 3 END ###


### Test full file ###
# if __name__ == "__main__":
#     pattern_path = "/Users/david/Sphere_files/Downsample replicate/output_csv/cluster_patterns.csv"
#     lc_path = "/Users/david/Sphere_files/Downsample replicate/output_csv/lowest common read depths - neuron - chr19.csv"
#     A = "/Users/david/Sphere_files/Downsample replicate/Clubcpg_re_run_July2021/sex_p35_neuron/male_p35_neuron.bam.chr19_cluster_results.csv"
#     B = "/Users/david/Sphere_files/Downsample replicate/Clubcpg_re_run_July2021/sex_p12_neuron/male_p12_neuron.bam.chr19_cluster_results.csv"
#     club_idx, clubcpg = downsampe_prep(pattern_path, lc_path, A, B)
#     print("finish prep")
#     print(clubcpg.shape)
#     # print("clubcpg",clubcpg)
#     output_df = run_downsample(clubcpg, club_idx)
#     print("finish running")
#     print("output", output_df)
#     output_df.to_csv("/Users/david/" + "output.csv")
#     print("time:", time.time() - start_time, "s")

### Test full file end ###
#
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("pattern_path", help="absolute path to new class labels", default=None)
    parser.add_argument("lc_path", help="absolute path to minimum read depth at each bin ", default=None)
    parser.add_argument("-A", help="absolute path to the first clubcpg cluster output", default=None)
    parser.add_argument("-B", help="absolute path to the second clubcpg cluster output", default=None)
    parser.add_argument("-sampleA", "--sampleA",
                        help=" when using 2 file mode, specify which sample in cluster output A to use, allowed inputs are A or B", default = None)
    parser.add_argument("-sampleB", "--sampleB",
                        help=" when using 2 file mode, specify which sample in cluster output A to use, allowed inputs are A or B", default=None)
    parser.add_argument("-chr", "--chromosome", help="Optional, perform only on one chromosome. ",default=None)
    parser.add_argument("-ncore","--ncore", help="the number of cores the downsampling can use", default=None)
    parser.add_argument("-o", "--output", help="folder to save imputed coverage data", default=None)
    parser.add_argument("-name", "--name", help="desired output file name", default="/output1.csv")
    args = parser.parse_args()
    # TODO: the -chr argument is not the used yet, think of how to use it.
    # Set output dir
    ncore = args.ncore
    if not args.output:
        output_folder = os.path.dirname(args.lc_path)
    else:
        output_folder = args.output
    try:
        os.mkdir(output_folder)
    except FileExistsError:
        print("Output folder already exists... no need to create it...")

    if not args.B: # if single clubcpg output:
        club_idx, clubcpg = downsampe_prep_single_file(args.pattern_path, args.lc_path, args.A)
        output_df = run_downsample(clubcpg,club_idx)
        output_df.to_csv(output_folder + args.name)
        print("time:", time.time() - start_time, "s")
    else:
        club_idx, clubcpg = downsampel_prep(args.pattern_path, args.lc_path, args.A, args.B, args.sampleA,args.sampleB)
        # TODO: test prep
        output_df = run_downsample(clubcpg, club_idx)
        # os.mkdir(output_folder + str(date.today()))
        output_df.to_csv(output_folder + args.name)
        # print(output_df)
        print("time:", time.time() - start_time, "s")
