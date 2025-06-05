import numpy as np
import random
import os
import time
from satlp import HybridSolver, SATasLPFeasibility

np.random.seed(239)

if __name__ == "__main__":
    n_vars = 20
    cnf_idx = 10
    original_f = f"cnfs/uf{n_vars}/uf{n_vars}-0{cnf_idx}.cnf"
    filename = f"cnfs/uf{n_vars}/uf{n_vars}-0{cnf_idx}.cnf.copy"
    os.system(f"cp {original_f} {filename}")
    solutions = f"solutions/uf{n_vars}/uf{n_vars}-0{cnf_idx}.cnf.out"
    # original_f = f"cnfs/cnf_sat.cnf"
    # filename = f"cnfs/cnf_sat.cnf.copy"
    # os.system(f"cp {original_f} {filename}")
    # solutions = f"solutions/cnf_sat.out"

    with open(solutions) as f:
        sols = [l for l in f.readlines()]

    hyb_solver = HybridSolver(filename, SATasLPFeasibility)
    witness = hyb_solver.solve()
    hyb_solver.verify(witness)

    sat_solver = BooleanSolver(original_f, verbose=0)
    sat_solver.solve()

