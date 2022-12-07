import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import statistics

p_val = [0.939043574058575, 0.6197887582883937,\
         0.30685953932470383, 0.30685953932470383,\
         0.939043574058575, 0.6197887582883937,\
         0.939043574058575, 0.939043574058575,
         0.6197887582883937, -0.0, 0.6197887582883937,
         0.30685953932470383, 0.6197887582883937,
         0.939043574058575, 0.30685953932470383,
         -0.0, 0.939043574058575, 0.939043574058575,
         0.30685953932470383, 0.6197887582883937,
         0.939043574058575, 0.30685953932470383,
         0.30685953932470383, 0.6197887582883937,
         -0.0, 0.6197887582883937, 0.939043574058575,
         0.6197887582883937, 0.6197887582883937,
         0.939043574058575, 0.6197887582883937,
         0.6197887582883937, 0.6197887582883937,
         -0.0, -0.0, 0.939043574058575, 0.6197887582883937, 0.6197887582883937, 0.6197887582883937, 0.6197887582883937, 0.30685953932470383, 0.30685953932470383, -0.0, 0.30685953932470383, 0.6197887582883937, -0.0, 0.6197887582883937, 0.6197887582883937, 0.6197887582883937, -0.0, 0.6197887582883937, -0.0, 0.939043574058575, 0.6197887582883937, 0.30685953932470383, 0.6197887582883937, 0.6197887582883937, 0.6197887582883937, 0.6197887582883937, 0.6197887582883937, 0.6197887582883937, 0.939043574058575, 0.939043574058575, 0.6197887582883937, 0.6197887582883937, 0.30685953932470383, 0.30685953932470383, 0.30685953932470383, 0.30685953932470383, 0.6197887582883937, 0.939043574058575, 0.30685953932470383, 0.30685953932470383, 0.30685953932470383, 0.6197887582883937, 0.30685953932470383, 0.6197887582883937, 0.30685953932470383, 0.30685953932470383, 0.30685953932470383, 0.6197887582883937, -0.0, -0.0, 0.6197887582883937, 0.6197887582883937, 0.30685953932470383, 0.6197887582883937, 0.939043574058575, 0.6197887582883937, 0.6197887582883937, 0.30685953932470383, 0.30685953932470383, 0.30685953932470383, 0.6197887582883937, 0.939043574058575, 0.6197887582883937, 0.6197887582883937, 0.6197887582883937, 0.939043574058575, 0.30685953932470383]

p_val = [1.9959527597159825, 1.7183014331322182, 0.9213498354804726, 1.7183014331322182, 1.7183014331322182, 1.4461163290995296,
 1.1801316432262405, 0.9213498354804726, 1.9959527597159825, 1.7183014331322182, 1.4461163290995296, 0.6711850122116552,
 1.4461163290995296, 0.9213498354804726, 1.4461163290995296, 1.9959527597159825, 2.278518516821683, 1.7183014331322182,
 1.4461163290995296, 1.4461163290995296, 1.7183014331322182, 1.4461163290995296, 2.8567929263111465, 1.9959527597159825,
 2.5655753326769597, 2.278518516821683, 0.9213498354804726, 1.9959527597159825, 2.278518516821683, 2.5655753326769597,
 2.8567929263111465, 1.9959527597159825, 1.4461163290995296, 1.4461163290995296, 2.278518516821683, 1.4461163290995296,
 1.7183014331322182, 1.4461163290995296, 1.1801316432262405, 1.7183014331322182, 1.1801316432262405, 1.1801316432262405,
 2.5655753326769597, 2.278518516821683, 1.9959527597159825, 1.9959527597159825, 1.1801316432262405, 1.4461163290995296,
 1.4461163290995296, 1.9959527597159825, 1.7183014331322182, 2.5655753326769597, 1.9959527597159825, 1.7183014331322182,
 1.7183014331322182, 1.9959527597159825, 1.9959527597159825, 1.7183014331322182, 1.4461163290995296, 1.4461163290995296,
 1.1801316432262405, 1.7183014331322182, 1.7183014331322182, 1.9959527597159825, 2.8567929263111465, 2.278518516821683,
 0.9213498354804726, 1.4461163290995296, 1.1801316432262405, 3.151910196397295, 1.4461163290995296, 1.1801316432262405,
 1.4461163290995296, 1.4461163290995296, 1.7183014331322182, 1.9959527597159825, 0.6711850122116552, 1.9959527597159825,
 1.4461163290995296, 1.7183014331322182, 1.4461163290995296, 1.9959527597159825, 2.278518516821683, 1.9959527597159825,
 2.5655753326769597, 1.4461163290995296, 0.9213498354804726, 1.7183014331322182, 0.9213498354804726, 1.7183014331322182,
 2.278518516821683, 0.6711850122116552, 1.9959527597159825, 1.9959527597159825, 1.4461163290995296, 1.7183014331322182,
 1.9959527597159825, 1.1801316432262405, 1.1801316432262405, 1.4461163290995296]

bin_name = "chr1_899823989"
class_label = 12
p_series = pd.Series(p_val)


def make_plot(bin_name,class_label,p_vals):
    mean = statistics.mean(p_vals)
    std = statistics.stdev(p_vals)

    plt.hist(p_vals, bins=20)
    plt.axvline(mean,color="k", linestyle = "dashed",label='{0:.4f}'.format(mean))
    plt.axvline(mean+std,color="y", linestyle = "dashed",label='{0:.4f}'.format(mean+std))
    plt.axvline(mean-std,color="y", linestyle = "dashed",label='{0:.4f}'.format(mean-std))
    plt.xticks(np.arange(-0.5, 3.5, 0.5))
    plt.legend(loc='upper right')

    plt.gca().set(title="bin: "+bin_name+ " class_label: " +str(class_label), xlabel="p_val",ylabel='Frequency')
    fig = plt.gcf()
    return fig





