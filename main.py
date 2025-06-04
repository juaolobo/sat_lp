from cnf_loader import CNFLoader
from linear_solver import (
    SATasLPOptimization,
    SATasLPFeasibility,
)
from boolean_solver import Solver
import numpy as np
from scipy.optimize import linprog
import random
import os
import time

np.random.seed(239)

def check_solutions(solutions, solution):
    errors = []
    min_err = n_vars
    argmin = -1
    for l in solutions:
        real_solution = np.array([int(x) for x in l.split()[:-1]])
        error = sum(real_solution != solution) - sum(solution == 0)
        errors.append(error)
        if error < min_err:
            min_err = error
            min_sol = real_solution
            wrong_idxs = real_solution[(real_solution != solution) & (solution != 0)]

    return min_err, min_sol, wrong_idxs

n_vars = 20
cnf_idx = 10
original_f = f"../cnfs/uf{n_vars}/uf{n_vars}-0{cnf_idx}.cnf"
filename = f"../cnfs/uf{n_vars}/uf{n_vars}-0{cnf_idx}.cnf.copy"
os.system(f"cp {original_f} {filename}")
solutions = f"../solutions/uf{n_vars}/uf{n_vars}-0{cnf_idx}.cnf.out"

# original_f = f"cnfs/cnf_sat.cnf"
# filename = f"cnfs/cnf_sat.cnf.copy"
# os.system(f"cp {original_f} {filename}")
# solutions = f"solutions/cnf_sat.out"

with open(solutions) as f:
    sols = [l for l in f.readlines()]

fixing = {}
last_integral = None
integral_vars = None

witness = np.zeros(n_vars)
it = 0
start = time.time()
restarts = 0

lp = SATasLPFeasibility(fixing=fixing, filename=filename, method='highs-ipm')
solver = Solver(filename, verbose=0)

while 1:

    lp.create_lp()
    last_witness = witness
    witness = lp.solve()

    # i.e. INFEASIBLE
    if witness == None:
        witness = [0.5 for _ in range(lp.n_vars())]
        restarts += 1
        break

    it += 1

    if lp.verify(witness) == True:
        break

    n = lp.n_vars()
    solution = np.zeros(n)
    new_fixing = {}
    for i, xi in enumerate(witness):
        if xi.is_integer():
            new_fixing[i+1] = xi
            solution[i] = i+1 if xi > 0 else -i-1    

    linear_sol = [xi if new_fixing[xi] else -xi for xi in new_fixing.keys()]
    last_integral = integral_vars
    integral_vars = sum(solution != 0)

    min_err, min_sol, wrong_idxs = check_solutions(sols, solution)

    print(f"Program achieved {integral_vars} integral variables")
    print(f"Partial LP solution: {linear_sol}")
    print(f"Solution is wrong {min_err} by variables {wrong_idxs}")
    print("-------------------------------------------------------")

    if new_fixing == fixing:

        solver.restart()
        resolved, formula = solver.propagate_linear(linear_sol)
        new_clauses = [f.clause for f in formula.formula[lp.m_clauses():]]

        if len(new_clauses) == 0:
            resolved, formula = solver.extend_solution()
            new_clauses = [f.clause for f in formula.formula[lp.m_clauses():]]

        for c in new_clauses:
            lp.add_clause(c)
            solver.list_clause.append(c)

        fixing = {abs(xi): 1.0 if xi > 0 else 0.0 for i, xi in enumerate(resolved)}
                
    else:
        fixing = new_fixing

    lp.restart(fixing=fixing)
    
if witness:
    sol = lp.verify(witness)
    if sol:
        print("SATISFIABLE")
    partial = [xi for xi in witness if xi.is_integer()]
    print(f"Achieved {len(partial)} integral variables")
    solution = [i+1 if witness[i] > 0 else -i-1 for i in range(len(partial))]
    print(solution)

else:
    breakpoint()

stop = time.time()
print(f"Iterations executed: {it}; elapsed time: {stop-start:.3f} seconds")

solver = Solver(original_f, verbose=0)
solver.solve()

