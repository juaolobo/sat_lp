import pandas as pd
from ast import literal_eval
from math import comb
from matplotlib import pyplot as plt
import sys
import numpy as np

def addlabels(x,y, ax):
    delta = 2
    props = dict(boxstyle='round', facecolor='beige', alpha=0.5)
    for i in range(len(x)):
        ax.text(i, y[i] + delta, round(y[i]), ha = 'center', parse_math=True, bbox=props)


ns = np.array([20, 50, 100, 200])
ticks = np.arange(len(ns))

opt_type = ["feasibility", "optimization"]
offset = -1
fig, ax = plt.subplots()

width = 0.25
mult = 0
ax.set_xticks(ticks + width, ns)

for t in opt_type:
    means = []
    stds = []
    n_satisfieds = []
    for n in ns:
        filename = f"experiments/data/integral_conv/uf{n}-{t}.csv"
        figsize=(13, 8)
        df = pd.read_csv(filename)

        mean_int = df["int_vars"].mean()
        std_int = df["int_vars"].std()
        n_satisfied = df["satisfied"].sum()

        means.append(mean_int)
        stds.append(std_int)
        n_satisfieds.append(n_satisfied)

    offset = width * mult
    mult += 1
    ax.bar(ticks+offset, means, width, yerr=stds)

ax.set_ylim([0, 100])
ax.set_title(f"Number of integral variables achieved in linear program")
ax.legend(opt_type, loc='upper left')
plt.savefig(f"figs/n_integral_vars_ipm.png")