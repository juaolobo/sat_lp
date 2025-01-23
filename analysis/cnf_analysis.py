from satlp import CNFLoader
import pandas as pd
import numpy as np
import sys

if len(sys.argv) < 2:
    exit(0)

filename = f"cnfs/{sys.argv[1]}.cnf"
cnf_handler = CNFLoader(filename)

sr = pd.Series(data=cnf_handler.clauses)
count = sr.explode().value_counts()
count_pos = count[count.index > 0]
count_neg = count[count.index < 0]
count_neg.index = abs(count_neg.index)

df = pd.DataFrame(columns=["appearances", "appearances_pos", "appearances_neg", "neighbors", "a/n ratio"])

df["appearances_pos"] = count_pos
df["appearances_neg"] = count_neg
df["appearances"] = df["appearances_pos"] + df["appearances_neg"]

df = df.sort_values("appearances")
# df.index = count_pos.index


sorted_max_neighbors = {}
clauses = cnf_handler.clauses
neighbors = {x : set() for x in range(1, cnf_handler.n_vars+1)}
for c in clauses:
    for l in c:
        if abs(l):
            a = set([abs(x) for x in set(c) - set([l])])
            neighbors[abs(l)] = neighbors[abs(l)] | a

sr_n = pd.Series(neighbors).apply(len)
# idxmax = sr_n.idxmax()
# max_v = sr_n.max()

# sorted_max_neighbors[idxmax] = max_v
# clauses = [c for c in clauses if idxmax not in np.abs(c)]


df["neighbors"] = sr_n
df["a/n ratio"] = df["appearances"]/df["neighbors"]

df = df.sort_values(["neighbors"])
print(df)

breakpoint()