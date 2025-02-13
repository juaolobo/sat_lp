import pandas as pd
from ast import literal_eval
from math import comb
import sys
import os

if len(sys.argv) < 2:
    exit(0)

dir_path = f"./figs/{sys.argv[1]}"
if not os.path.exists(dir_path):
    os.makedirs(dir_path)

filename = f"experiments/data/{sys.argv[1]}.csv"
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
df_vars_ncvg.to_csv(f"experiments/processed/{sys.argv[1]}-non_convergent_vars.csv")

fix_len = fixed_vars.apply(len)
ranked_len = pd.DataFrame(fix_len.value_counts())

max_comb = max(ranked_len.index)
ranked_len["comb_total"] = [comb(max_comb, i) for i in ranked_len["count"].index]
ranked_len_pct = ranked_len["count"]/ranked_len["comb_total"]
ranked_len["pct"] = ranked_len_pct
df_len_ncvg = (ranked_len
                .sort_values(by=["pct", "fixing"])
                .drop(columns="comb_total"))

df_len_ncvg.to_csv(f"experiments/processed/{sys.argv[1]}-non_convergent_len.csv")





