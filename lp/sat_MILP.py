"""
What we are going to do is to implement the original idea for 3-SAT via LP,
but we the constraints will be set to the linear system of the Fourier expansions
of the disjunctions

x1 v x2 v x3 = max x1 x2 x3 = max x1 (max x2 x3) = 1/2 + x1/2 + (max x2 x3)/2 - x1 max(x2 x3)/2
= 1/2 + x1/2 + (1/2 + x2/2 + x3/2 - x2x3/2)/2 - x1 (1/2 + x2/2 + x3/2 - x2x3/2) / 2
= 1/2 + x1/2 + 1/4 + x2/4 + x3/4 - x2x3/4 - x1/4 - x1x2/4 - x1x3/4 + x1x2x4/4
x1 v x2 v x3 = 3/4 + x1/4 + x2/4 + x3/4 - (x2x3 + x1x2 + x1x3)/4 + x1x2x3/4 = 1
x1/4 + x2/4 + x3/4 - (x2x3 + x1x2 + x1x3)/4 + x1x2x3/4 = 1/4
"""
""" Read input from a file in DMACS format """

from ortools.linear_solver import pywraplp

import numpy as np
import pandas as pd
class Cnf_lp:

    def __init__(self):
        self.n_vars = None
        self.m_clauses = None
        self.clauses = None
        self.solver   = pywraplp.Solver.CreateSolver("SCIP")
        if not self.solver:
            raise Exception("Solver creation failed")

    def _load_from_file(self, filename):

        with open(filename) as f:
            cnf_data = [line.rstrip() for line in f if line != '\n']
            m_clauses = int(cnf_data[0].split()[-1])
            n_variables = int(cnf_data[0].split()[-2])
            clauses = [
                [int(i) for i in l.split()[:-1]]
                    for l in cnf_data[1:]
            ]

        self.n_vars = n_variables
        self.m_clauses = m_clauses
        self.clauses = clauses


    def _init_objects(self):

        M = 1
        self.vars = [self.solver.NumVar(0, 1.0, f"x_{i}") for i in range(1, self.n_vars + 1)]
        self.vars_prime = [self.solver.NumVar(0, 0.5, f"x_{i}") for i in range(1, self.n_vars + 1)]
        self.int_vars = [self.solver.BoolVar(f"b_{i}") for i in range(1, self.n_vars + 1)]

        res = [1 for c in self.clauses]

        for i in range(self.m_clauses):
            res[i] = np.sum([0 if x > 0 else -1 for x in self.clauses[i]]) + 1

        # create clause constraints
        for i, c in enumerate(self.clauses):
            coefs = []
            for l in c:
                l_v = np.abs(l)
                sgn = l_v/l
                idx = l_v - 1
                coefs.append(sgn*self.vars[idx])
            
            self.solver.Add(self.solver.Sum(coefs) >= res[i])

        # add absolute value contraints
        for i in range(self.n_vars):

            self.solver.Add(self.vars[i]-1/2 + M*self.int_vars[i] >= self.vars_prime[i])
            self.solver.Add(-self.vars[i]+1/2 + M*(1 - self.int_vars[i]) >= self.vars_prime[i])
            self.solver.Add(self.vars[i]-1/2 <= self.vars_prime[i])
            self.solver.Add(-self.vars[i]+1/2 <= self.vars_prime[i])

    def _create_optimization(self):
        # optimizing for OBJ(x) = 1, for now            
        self.solver.Maximize(self.solver.Sum(self.vars_prime))

    def create_lp(self, filename=None):
        self._load_from_file(filename)
        self._init_objects()
        self._create_optimization()


if __name__ == "__main__":

    lp_obj = Cnf_lp()
    # filename = "cnf_sat.txt"
    filename = "uf20-018.cnf"
    lp_obj.create_lp(filename)
    lp_obj.solver.Solve()
    witness = {v.name(): v.solution_value() for v in lp_obj.vars}
    print(witness)
    print("----------------")
    result = lp_obj.solver.Objective().Value()
    print(result)
