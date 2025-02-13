import pandas as pd
from ast import literal_eval
from math import comb
from matplotlib import pyplot as plt
import sys
import os

def addlabels(x,y):
    delta = 2e-2
    props = dict(boxstyle='round', facecolor='beige', alpha=0.5)
    for i in range(len(x)):
        plt.text(i, y[i] + delta, round(y[i], 2), ha = 'center', parse_math=True, bbox=props)

if len(sys.argv) < 2:
    exit(0)

dir_path = f"./figs/{sys.argv[1]}"
if not os.path.exists(dir_path):
    os.makedirs(dir_path)

filename_vars = f"experiments/processed/{sys.argv[1]}-non_convergent_vars.csv"
filename_len = f"experiments/processed/{sys.argv[1]}-non_convergent_len.csv"
figsize=(13, 8)
df_vars = pd.read_csv(filename_vars).set_index("fixing")
df_len = pd.read_csv(filename_len).set_index("n_fixed")

fig = plt.figure(figsize=figsize)
df_vars["pct"].plot(kind='bar', title='Percentage of times that fixing variables lead to divergent solution')
addlabels(df_vars.index, df_vars["pct"].values)
plt.savefig(f"{dir_path}/{sys.argv[1]}-non_convergent_vars.png")

fig = plt.figure(figsize=figsize)
df_len["pct"].plot(kind='bar', title="Number of fixed variables in divergent solution")
addlabels(df_len.index, df_len["pct"].values)
plt.savefig(f"{dir_path}/{sys.argv[1]}-non_convergent_len.png")


