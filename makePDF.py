import itertools
import math
import statistics
from functools import partial

# import dask.dataframe as dd
import pandas as pd
import numpy as np
import time
import scipy as sp
from scipy import stats
from modify_cluster_output_file import *
from multiprocessing import Pool, cpu_count
# from dask.distributed import LocalCluster, Client
import multiprocessing
from typing import Callable, Tuple, Union
from datetime import date
import matplotlib.pyplot as plt
from make_plot import make_plot
from matplotlib.backends.backend_pdf import PdfPages


with PdfPages("/Users/david/Sphere_files/Downsample replicate/CluBCpG demos/output_csv/outplot3.pdf") as pdf:
    for i in range(10):
        p_vals = [1.9959527597159825, 1.7183014331322182, 0.9213498354804726, 1.7183014331322182,
                  1.7183014331322182, 1.4461163290995296,
                  1.1801316432262405, 0.9213498354804726, 1.9959527597159825, 1.7183014331322182,
                  1.4461163290995296, 0.6711850122116552,
                  1.4461163290995296, 0.9213498354804726, 1.4461163290995296, 1.9959527597159825,
                  2.278518516821683, 1.7183014331322182,
                  1.4461163290995296, 1.4461163290995296, 1.7183014331322182, 1.4461163290995296,
                  2.8567929263111465, 1.9959527597159825,
                  2.5655753326769597, 2.278518516821683, 0.9213498354804726, 1.9959527597159825,
                  2.278518516821683, 2.5655753326769597,
                  2.8567929263111465, 1.9959527597159825, 1.4461163290995296, 1.4461163290995296,
                  2.278518516821683, 1.4461163290995296,
                  1.7183014331322182, 1.4461163290995296, 1.1801316432262405, 1.7183014331322182,
                  1.1801316432262405, 1.1801316432262405,
                  2.5655753326769597, 2.278518516821683, 1.9959527597159825, 1.9959527597159825,
                  1.1801316432262405, 1.4461163290995296,
                  1.4461163290995296, 1.9959527597159825, 1.7183014331322182, 2.5655753326769597,
                  1.9959527597159825, 1.7183014331322182,
                  1.7183014331322182, 1.9959527597159825, 1.9959527597159825, 1.7183014331322182,
                  1.4461163290995296, 1.4461163290995296,
                  1.1801316432262405, 1.7183014331322182, 1.7183014331322182, 1.9959527597159825,
                  2.8567929263111465, 2.278518516821683,
                  0.9213498354804726, 1.4461163290995296, 1.1801316432262405, 3.151910196397295,
                  1.4461163290995296, 1.1801316432262405,
                  1.4461163290995296, 1.4461163290995296, 1.7183014331322182, 1.9959527597159825,
                  0.6711850122116552, 1.9959527597159825,
                  1.4461163290995296, 1.7183014331322182, 1.4461163290995296, 1.9959527597159825,
                  2.278518516821683, 1.9959527597159825,
                  2.5655753326769597, 1.4461163290995296, 0.9213498354804726, 1.7183014331322182,
                  0.9213498354804726, 1.7183014331322182,
                  2.278518516821683, 0.6711850122116552, 1.9959527597159825, 1.9959527597159825,
                  1.4461163290995296, 1.7183014331322182,
                  1.9959527597159825, 1.1801316432262405, 1.1801316432262405, 1.4461163290995296]
        mean = statistics.mean(p_vals)
        std = statistics.stdev(p_vals)

        fig = plt.figure(figsize=(10, 10))
        plt.hist(p_vals, bins=20)
        plt.axvline(mean, color="k", linestyle="dashed", label='{0:.4f}'.format(mean))
        plt.axvline(mean + std, color="y", linestyle="dashed", label='{0:.4f}'.format(mean + std))
        plt.axvline(mean - std, color="y", linestyle="dashed", label='{0:.4f}'.format(mean - std))
        plt.xticks(np.arange(-0.5, 3.5, 0.5))
        plt.legend(loc='upper right')

        plt.gca().set(title="bin: " + "chr1" + " class_label: " + str(5), xlabel="p_val",
                      ylabel='Frequency')

        # for key, item in grouped_df:
        #     global count
        #     count += 1
        #     if count <= 30:
        #         print(grouped_df.get_group(key)["p_val"])
        #         # draw yor plot
        #         # fig = make_plot(bin_name,class_label=key, p_vals=grouped_df.get_group(key)["p_val"].tolist())
        #
        #         p_vals = grouped_df.get_group(key)["p_val"].tolist()
        #         mean_check = statistics.mean(p_vals)
        #         mean = grouped_df.get_group(key)["p_val"].mean()
        #         std = grouped_df.get_group(key)["p_val"].median()
        #         median = grouped_df.get_group(key)["p_val"].std()
        #
        #         # fig = plt.figure(figsize=(10, 10))
        #         plt.hist(p_vals, bins=20)
        #         plt.axvline(mean, color="k", linestyle="dashed", label='{0:.4f}'.format(mean))
        #         plt.axvline(mean + std, color="y", linestyle="dashed", label='{0:.4f}'.format(mean + std))
        #         plt.axvline(mean - std, color="y", linestyle="dashed", label='{0:.4f}'.format(mean - std))
        #         plt.axvline(median, color="r", linestyle="dashed", label='{0:.4f}'.format(median))
        #         plt.xticks(np.arange(-0.5, 3.5, 0.5))
        #         plt.legend(loc='upper right')
        #
        #         plt.gca().set(title="bin: " + bin_name + " class_label: " + str(key), xlabel="p_val",
        #                       ylabel='Frequency')
        #         # figs.append(fig)
        #         # plt.close()
        pdf.savefig(fig)