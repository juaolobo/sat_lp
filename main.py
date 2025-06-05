import numpy as np
import random
import os
import time
from satlp import HybridSolver, SATasLPFeasibility, BooleanSolver

np.random.seed(239)

if __name__ == "__main__":
    n_vars = 20
    cnf_idx = 10
    filename = f"formulas/cnfs/uf{n_vars}/uf{n_vars}-0{cnf_idx}.cnf"

    hyb_solver = HybridSolver(filename, SATasLPFeasibility)
    witness = hyb_solver.solve()
    hyb_solver.verify(witness)

    sat_solver = BooleanSolver(filename, verbose=0)
    sat_solver.solve()

