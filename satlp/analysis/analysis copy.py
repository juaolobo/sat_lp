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

filename = f"experiments/{sys.argv[1]}-experiments.csv"
figsize=(13, 8)
df = pd.read_csv(filename)
cvg = df[df["is_solution"] == True]
ncvg = df[df["is_solution"] == False]

# rank vars that lead to most trouble after relaxation
fixed_vars = ncvg["fixing"].apply(literal_eval)
fixed_vars = fixed_vars.apply(lambda x: x.keys())
ranked_vars = pd.DataFrame(fixed_vars.explode().value_counts())
ranked_vars_pct = ranked_vars/len(fixed_vars)
ranked_vars["pct"] = ranked_vars_pct

df_vars_ncvg = ranked_vars

fig = plt.figure(figsize=figsize)
df_vars_ncvg["pct"].plot(kind='bar', title='Percentage of times that variables appeared in divergent solution')
addlabels(df_vars_ncvg["pct"].index, df_vars_ncvg["pct"].values)
plt.savefig(f"figs/{sys.argv[1]}-non_convergent_vars.png")

fix_len = fixed_vars.apply(len)
ranked_len = pd.DataFrame(fix_len.value_counts())

max_comb = max(ranked_len.index)
ranked_len["comb_total"] = [comb(max_comb, i) for i in ranked_len["count"].index]
ranked_len_pct = ranked_len["count"]/ranked_len["comb_total"]
ranked_len["pct"] = ranked_len_pct

breakpoint()
df_len_ncvg = (ranked_len
                .sort_values(by=["pct", "fixed_vars"])
                .drop(columns="comb_total"))

fig = plt.figure(figsize=figsize)
df_len_ncvg["pct"].plot(kind='bar', title="Number of relaxed variables in divergent solution")
addlabels(df_len_ncvg["pct"].index, df_len_ncvg["pct"].values)
plt.savefig(f"figs/{sys.argv[1]}-non_convergent_len.png")

breakpoint()





