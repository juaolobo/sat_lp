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

    def __init__(self, fix_variables = {}):
        self.n_vars = None
        self.m_clauses = None
        self.clauses = None
        self.solver   = pywraplp.Solver.CreateSolver("GLOP")
        self.fix_variables = fix_variables
        if not self.solver:
            raise Exception("Solver creation failed")

    def _load_from_file(self, filename):

        with open(filename) as f:
            cnf_data = f.readlines()
            m_clauses = int(cnf_data[0].split()[-1])
            n_variables = int(cnf_data[0].split()[-2])
            clauses = [[int(i) for i in l.split()[:-1]] for l in cnf_data[1:]]

        self.n_vars = n_variables
        self.m_clauses = m_clauses
        self.clauses = clauses
        if self.fix_variables:
            self.free_vars = list(set(self.fix_variables.keys()) - set(range(self.n_vars)))
        else:
            self.free_vars = range(self.n_vars)


    def _init_objects(self):

        self.vars = [self.solver.NumVar(0, 1.0, f"x_{i}") for i in range(1, self.n_vars + 1)]
        self.pos_vars = [self.solver.NumVar(0, 0.5, f"x_{i}_p") for i in range(1, self.n_vars + 1)]
        self.neg_vars = [self.solver.NumVar(0, 0.5, f"x_{i}_n") for i in range(1, self.n_vars + 1)]

        # for idx, v in self.fix_variables.items():
        #     self.vars[idx].SetBounds(v, v)
        #     self.pos_vars[idx].SetBounds(0,0)
        #     self.neg_vars[idx].SetBounds(0,0)

        res = [1 for c in self.clauses]

        for i in range(self.m_clauses):
            res[i] = np.sum([0 if x > 0 else -1 for x in self.clauses[i]]) + 1

        constrs = [self.solver.Constraint(res[i].item(), self.solver.infinity()) for i in range(self.m_clauses)]

        for i, c in enumerate(self.clauses):
            for l in c:
                l_v = np.abs(l)
                sgn = l_v/l
                idx = l_v - 1
                constrs[i].SetCoefficient(self.vars[idx], sgn)

        for i in range(self.n_vars):

            # x_p + x_n = x - 1/2 => x_p + x_n -x = -1/2 => x - x_p - x_n = 1/2
            # this does not work because of maximization :C
            # try mixed integer programming and relaxations
            # constr = self.solver.Constraint(1/2, 1/2)
            # constr.SetCoefficient(self.vars[i], 1)
            # constr.SetCoefficient(self.pos_vars[i], -1)
            # constr.SetCoefficient(self.neg_vars[i], 1)
            # constrs.append(constr)
        
            # constr = self.solver.Constraint(0.0, 0.5)
            # constr.SetCoefficient(self.pos_vars[i], 1)
            # constr.SetCoefficient(self.neg_vars[i], 1)
            # constrs.append(constr)

    def _create_optimization(self):
        # optimizing for OBJ(x) = 1, for now
        objective = self.solver.Objective()
        for i in range(self.n_vars):
            # objective.SetCoefficient(self.pos_vars[i], 1)
            # objective.SetCoefficient(self.neg_vars[i], -1)

        objective.SetMaximization()

    def create_lp(self, filename=None):
        self._load_from_file(filename)
        self._init_objects()
        self._create_optimization()

    def extract_fix_variables(self):

        fix_variables = {}
        for idx, v in enumerate(self.vars):
            if v.solution_value() == 1 or v.solution_value() == 0:
                fix_variables[idx] = v.solution_value()

        return fix_variables


if __name__ == "__main__":

    fix_variables = {}
    lp_obj = Cnf_lp(fix_variables = fix_variables)
    filename = "cnf_sat.txt"
    # filename = "uf20-018.cnf"
    lp_obj.create_lp(filename)
    lp_obj.solver.Solve()
    witness = {v.name(): v.solution_value() for v in lp_obj.vars}
    print(witness)
    result = lp_obj.solver.Objective().Value()
    print(result)

    # for i in range(4):
    #     lp_obj = Cnf_lp(fix_variables = fix_variables)
    #     lp_obj.create_lp("uf20-018.cnf")
    #     lp_obj.solver.Solve()
    #     fix_variables = lp_obj.extract_fix_variables()
    #     str_fix = [f"x_{i+1} = {v}" for i,v in fix_variables.items()]
    #     witness = {v.name(): v.solution_value() for v in lp_obj.vars}
    #     print(witness)
    #     print(str_fix)
    #     print("----------------")
    # result = lp_obj.solver.Objective().Value()

    # witness_pos = {v.name(): v.solution_value() for v in lp_obj.pos_vars}
    # witness_neg = {v.name(): v.solution_value() for v in lp_obj.neg_vars}
    # print(witness)
    # print(witness_pos)
    # print(witness_neg)
    # print(result)