from sat_MILP import SATasLP
from tqdm import tqdm
from itertools import combinations
import csv
import sys

with open(f"experiments.txt", "w") as f:
    writer = csv.writer(f, delimiter=',', quoting=csv.QUOTE_ALL)
    columns = ['result','is_solution','witness','betas','relaxed_vars']
    writer.writerow(columns)
    for i in range(1, 21):
        cmbs = [c for c in combinations(range(1, 21), i)]
        print(f"Testing combinations for (20 {i})")
        for j in tqdm(range(len(cmbs))):
            cmb = cmbs[j]
            lp_obj = SATasLP(relaxed_vars=cmb)
            filename = "cnfs/uf20-018.cnf"
            lp_obj.create_lp(filename)
            res, witness, betas = lp_obj.solve()
            ok = lp_obj.verify(witness)
            csv_values = [res, ok, witness, betas, cmb]
            writer.writerow(csv_values)

            