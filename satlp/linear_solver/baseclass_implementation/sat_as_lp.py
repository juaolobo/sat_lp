from satlp.linear_solver.sat_baseclass.sat_baseclass import SATasLPBaseclass

from scipy.optimize import linprog
import numpy as np
from abc import ABC, abstractmethod

class SATasLP(SATasLPBaseclass):

    def __init__(self, filename=None, cnf_handler=None, method=None):
        super().__init__(filename, cnf_handler)
        self.solver = linprog
        self.g = lambda x: 1 - self.fixing[abs(x)] if x < 0 else self.fixing[abs(x)]
        self.A_ub = None
        self.y_ub = None
        self.A_eq = None
        self.y_eq = None
        self.c = None
        self.bounds = None
        self.method = method

        if not self.solver:
            raise Exception("Solver creation failed")

    def solve(self, x0=None):

        result = self.solver(
            self.c, 
            x0=x0,
            A_eq=self.A_eq, 
            b_eq=self.y_eq,
            A_ub=self.A_ub, 
            b_ub=self.y_ub, 
            bounds=self.bounds, 
            method=self.method
        )
        x = result.x

        if result.success:
            witness = [x[i-1].item() if i not in self.fixing.keys() else self.fixing[i] for i in range(1,len(x)+1)]
            return witness

        print("INFEASIBLE")

class SATasMILP(SATasLPBaseclass):

    def __init__(self, filename=None):
        super().__init__(filename)
        self.solver = pywraplp.Solver.CreateSolver("SCIP")
        if not self.solver:
            raise Exception("Solver creation failed")

    def solve(self):
        s = self.solver.Solve()

        if s == self.solver.INFEASIBLE:
            print("INFEASIBLE")
            self.solver.Clear()
            return s, [], []

        result = self.solver.Objective().Value()
        witness = self.round([v.solution_value() if not isinstance(v, int) else v for v in self.vars])

        return s, result, witness

        