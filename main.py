from satlp import CNFLoader
from satlp import (
    SATasLPOptimization,
    SATasLPFeasibility,
    SATasMILPOptimization,
    SATasMILPFeasibility
)
from cdcl import CDCL_Solver
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
cnf_idx = 4
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

# fixing = {1:0, 4: 0, 5: 1, 6: 1, 8: 0, 9: 0, 10: 1, 13: 0, 15: 0, 16: 1, 17: 1, 18: 1}
# ass = np.array([-1, -4, 5, 6, -8, -9, 10, -13, -15, 16, 17, 18])
# fixing = {1: 0, 4: 0, 5: 1, 6: 1, 8: 0, 9: 0, 10: 1, 13: 0, 15: 0, 16: 1, 17: 1, 18: 1}
# ass = np.array([1, 2, -3, 6, 11, -14, 16, 17, -19])
# fixing = {1: 1, 2: 1, 3: 0, 6: 1, 11: 1, 14: 0, 16: 1, 17: 1, 19:0}
fixing = {}
last_integral = None
integral_vars = None

witness = np.zeros(n_vars)
it = 0
start = time.time()
restarts = 0

lp = SATasLPFeasibility(fixing=fixing, filename=filename, method='highs-ipm')
solver = CDCL_Solver(filename, verbose=0)
solver.solve()

while 1:

    lp.create_lp()
    last_witness = witness
    witness = lp.solve()
    if witness == None:
        fixing = {}
        witness = [0.5 for _ in range(lp.n_vars())]
        restarts += 1
        break

    it += 1

    if lp.verify(witness) == True:
        break

    n = lp.n_vars()
    solution = np.zeros(n)
    new_fixing = {i+1: xi for i, xi in enumerate(witness) if xi.is_integer()}
    dmacs_sol = [xi if new_fixing[xi] else -xi for xi in new_fixing.keys()]

    for i, xi in enumerate(witness):
        if xi.is_integer():
            new_fixing[i+1] = xi
            solution[i] = i+1 if xi > 0 else -i-1

    last_integral = integral_vars
    integral_vars = sum(solution != 0)

    min_err, min_sol, wrong_idxs = check_solutions(sols, solution)

    print(f"Program achieved {integral_vars} integral variables")
    print(f"Partial LP solution: {dmacs_sol}")
    print(f"Solution is wrong {min_err} by variables {wrong_idxs}")
    print("-------------------------------------------------------")

    # if new_fixing == fixing:
    solver.restart(hypotheses=dmacs_sol)
    resolved, formula = solver.solve_for_real()
    new_clauses = [f.clause for f in formula.formula[lp.m_clauses():]]
    if len(new_clauses) > 0:
        for c in new_clauses:
            lp.add_clause(c)
            solver.list_clause.append(c)

        fixing = {abs(xi): 1 if xi > 0 else 0 for i, xi in enumerate(resolved)}
        # else:
        #     fixing = {}
            
    else:
        fixing = new_fixing

    breakpoint()    
    print(fixing)
    # if new_fixing == fixing and 0:
    #     # while dmacs_sol is not None:
    #     for _ in range(n*len(dmacs_sol)):
    #     # while resolved is None:
    #         solver.restart(hypotheses=dmacs_sol)
    #         resolved, formula = solver.solve_for_conflict2()
    #         if resolved is not None:
    #             print(f"Resolved solution: {resolved}")
    #             new_clauses = [f.clause for f in formula.formula[lp.m_clauses():]]
    #             print(f"Learned clauses:")
    #             for c in new_clauses:
    #                 print("\t" + str(c))
                    
    #             if len(new_clauses) > 0:
    #                 for c in new_clauses:
    #                     lp.add_clause(c)
    #                     solver.list_clause.append(c)

    #                 if len(resolved) > 0:
    #                     fixing = {abs(xi): 1 if xi > 0 else 0 for xi in resolved}
    #                 else:
    #                     fixing = {}

    #                 print(f"Total clauses: {lp.m_clauses()}")
    #                 print(f"Learned {len(new_clauses)} clause(s)")
    #                 print("-------------------------------------------------------")
    #         else:
    #             break

    #         # dmacs_sol=resolved

    # else:
    #     fixing = new_fixing

    lp.restart(fixing=fixing)
    
solver = CDCL_Solver(original_f, verbose=0, hypotheses=[])
solver.solve()

if witness:
    sol = lp.verify(witness)
    if sol:
        print("SATISFIABLE")
    partial = [xi for xi in witness if xi.is_integer()]
    print(f"Achieved {len(partial)} integral variables")
    solution = [i+1 if witness[i] > 0 else -i-1 for i in range(len(partial))]
    print(solution)
    # print(f"Restarts: {restarts}")

else:
    breakpoint()

stop = time.time()
print(f"Iterations executed: {it}; elapsed time: {stop-start:.3f} seconds")

