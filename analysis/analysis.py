import pandas as pd
from ast import literal_eval
from math import comb
from matplotlib import pyplot as plt
import sys

def addlabels(x,y):
    delta = 2e-2
    props = dict(boxstyle='round', facecolor='beige', alpha=0.5)
    for i in range(len(x)):
        plt.text(i, y[i] + delta, round(y[i], 2), ha = 'center', parse_math=True, bbox=props)

if len(sys.argv) < 2:
    exit(0)

filename = f"experiments/data/{sys.argv[1]}.csv"
figsize=(13, 8)
df = pd.read_csv(filename)
cvg = df[df["result"] == True]
ncvg = df[df["result"] == False]

breakpoint()
# rank vars that lead to most trouble after relaxation
relaxed_vars = ncvg["relaxed_vars"].apply(literal_eval)
ranked_vars = pd.DataFrame(relaxed_vars.explode().value_counts())
ranked_vars_pct = ranked_vars/len(relaxed_vars)
ranked_vars["pct"] = ranked_vars_pct

df_vars_ncvg = ranked_vars

fig = plt.figure(figsize=figsize)
df_vars_ncvg["pct"].plot(kind='bar', title='Percentage of times that variables appeared in divergent solution')
addlabels(df_vars_ncvg["pct"].index, df_vars_ncvg["pct"].values)
plt.savefig(f"figs/{sys.argv[1]}-non_convergent_vars.png")

relax_len = relaxed_vars.apply(len)
ranked_len = pd.DataFrame(relax_len.value_counts())

max_comb = max(ranked_len.index)
ranked_len["comb_total"] = [comb(max_comb, i) for i in ranked_len["count"].index]
ranked_len_pct = ranked_len["count"]/ranked_len["comb_total"]
ranked_len["pct"] = ranked_len_pct

df_len_ncvg = (ranked_len
                .sort_values(by=["pct", "relaxed_vars"])
                .drop(columns="comb_total"))

fig = plt.figure(figsize=figsize)
df_len_ncvg["pct"].plot(kind='bar', title="Number of relaxed variables in divergent solution")
addlabels(df_len_ncvg["pct"].index, df_len_ncvg["pct"].values)
plt.savefig(f"figs/{sys.argv[1]}-non_convergent_len.png")

breakpoint()





