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
from itertools import combinations

import sys

import numpy as np
import pandas as pd
class Cnf_lp:

    def __init__(self):
        self.n_vars = None
        self.m_clauses = None
        self.clauses = None
        self.solver   = pywraplp.Solver.CreateSolver("GLOP")
        if not self.solver:
            raise Exception("Solver creation failed")

    def _load_from_file(self, filename):

        with open(filename) as f:

            cnf_data = []
            for line in f:
                stripped = line.rstrip()
                if len(stripped) > 0:
                    cnf_data.append(stripped)

            m_clauses = int(cnf_data[0].split()[-1])
            n_variables = int(cnf_data[0].split()[-2])
            clauses = [[int(i) for i in l.split()[:-1]] for l in cnf_data[1:]]

        self.n_vars = n_variables
        self.m_clauses = m_clauses
        self.clauses = [sorted(c, key=abs) for c in clauses]


    def _init_objects(self):

        self.vars = [self.solver.NumVar(-1.0, 1.0, f"x_{i}") for i in range(1, self.n_vars + 1)]

        self.multivars2 = {}
        self.multivars3 = {}

        for c in self.clauses:
            constr = self.solver.Constraint(1/4, 1/4)

            # set coefficients for x1/4, x2/4, x3/4
            for v in c:
                idx = abs(v)-1
                constr.SetCoefficient(self.vars[idx], np.sign(v)/4)

            # we combinate based on abs value to make sure indexing is consistent
            idx_c = np.abs(c)
            idx_comb2 = [x for x in combinations(idx_c, 2)]
            for cmb in idx_comb2:
                self.multivars2[cmb] = self.solver.NumVar(-1.0, 1.0, f"x_{cmb[0]}.x_{cmb[1]}")

            # set coefficients for -x1x2/4, -x1x3/4, -x2x3/4
            comb2 = [x for x in combinations(c, 2)]
            for cmb2 in comb2:
                coef = np.sign(np.prod(cmb2))
                _cmb = [np.abs(x).item() for x in cmb2]
                cmb = tuple(_cmb)
                constr.SetCoefficient(self.multivars2[cmb], -coef/4)

            # set coefficient for +x1x2x3/4
            coef = np.sign(np.prod(c))
            _cmb = [x.item() for x in idx_c]
            cmb = tuple(_cmb)
            self.multivars3[cmb] = self.solver.NumVar(-1.0, 1.0, f"x_{cmb[0]}.x_{cmb[1]}.x_{cmb[2]}")
            constr.SetCoefficient(self.multivars3[cmb], coef/4)

            breakpoint()

    def _create_optimization(self):
        # optimizing for OBJ(x) = 1, for now
        obj = self.solver.Objective()
        obj.SetMaximization()


    def create_lp(self, filename=None):
        self._load_from_file(filename)
        self._init_objects()
        self._create_optimization()


if __name__ == "__main__":


    lp_obj = Cnf_lp()
    if len(sys.argv) > 1:
        lp_obj.create_lp(sys.argv[1])
    else:
        lp_obj.create_lp("cnf_sat.txt")

    lp_obj.solver.Solve()

    result = lp_obj.solver.Objective().Value()

    witness = {v.name(): v.solution_value() for v in lp_obj.vars}
    print(witness)
    print(result)