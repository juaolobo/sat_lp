from satlp.sat_baseclass.sat_baseclass import SATasLPBaseclass

from ortools.linear_solver import pywraplp
from scipy.optimize import linprog
import numpy as np
from abc import ABC, abstractmethod


class SATasLPSimplex(SATasLPBaseclass):

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

class SATasLPIP(SATasLPBaseclass):

    def __init__(self, filename=None):
        super().__init__(filename)
        self.solver = linprog
        self.A = None
        self.y = None
        self.c = None
        self.bounds = None
        self.method = None
        if not self.solver:
            raise Exception("Solver creation failed")

    def solve(self):

        return self.solver(self.c, A_ub=self.A, b_ub=self.y, bounds=self.bounds, method=self.method)