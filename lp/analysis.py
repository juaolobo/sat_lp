import pandas as pd
from ast import literal_eval
from math import comb
from matplotlib import pyplot as plt
import sys

if len(sys.argv) < 2:
    exit(0)

filename = f"experiments/{sys.argv[1]}-experiments.csv"

df = pd.read_csv(filename)
cvg = df[df["is_solution"] == True]
ncvg = df[df["is_solution"] == False]

# rank vars that lead to most trouble after relaxation
relaxed_vars = ncvg["relaxed_vars"].apply(literal_eval)
ranked_vars = pd.DataFrame(relaxed_vars.explode().value_counts())
ranked_vars_pct = ranked_vars/len(relaxed_vars)
ranked_vars["pct"] = ranked_vars_pct

df_vars_ncvg = ranked_vars

fig = plt.figure()
df_vars_ncvg["pct"].plot(kind='bar')
plt.savefig(f"figs/{sys.argv[1]}non_convergent_vars.png")

relax_len = relaxed_vars.apply(len)
ranked_len = pd.DataFrame(relax_len.value_counts())

max_comb = max(ranked_len.index)
ranked_len["comb_total"] = [comb(max_comb, i) for i in ranked_len["count"].index]
ranked_len_pct = ranked_len["count"]/ranked_len["comb_total"]
ranked_len["pct"] = ranked_len_pct

df_len_ncvg = (ranked_len
                .sort_values(by=["pct", "relaxed_vars"])
                .drop(columns="comb_total"))

fig = plt.figure()
df_len_ncvg["pct"].plot(kind='bar', title="num_relaxed_vars")
plt.savefig(f"figs/{sys.argv[1]}non_convergent_len.png")
breakpoint()

witness = ncvg["witness"].apply(literal_eval)

results = ncvg["result"].value_counts()




