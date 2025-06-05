import numpy as np
import random
import os
import time
from satlp import HybridSolver, SATasLPFeasibility, BooleanSolver

if __name__ == "__main__":
    n_vars = 100
    cnf_idx = 4
    filename = f"formulas/cnfs/uf{n_vars}/uf{n_vars}-0{cnf_idx}.cnf"

    start = time.time()
    hyb_solver = HybridSolver(filename, SATasLPFeasibility)
    witness = hyb_solver.solve()
    stop = time.time()
    print(f"Elapsed time: {start - stop}s")
    hyb_solver.verify(witness)

    sat_solver = BooleanSolver(filename, verbose=0)
    sat_solver.solve()

