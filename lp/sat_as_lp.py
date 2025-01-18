from ortools.linear_solver import pywraplp
from cnf_loader import CNFLoader
import numpy as np
from abc import ABC, abstractmethod

class SATasLP(ABC):

    def __init__(self, filename=None, relaxed_vars=[]):
        self.cnf_handler = CNFLoader(filename)
        self.relaxed_vars = relaxed_vars
        self.solver = pywraplp.Solver.CreateSolver("SCIP")
        if not self.solver:
            raise Exception("Solver creation failed")

    @abstractmethod
    def _init_objects(self):
        pass

    @abstractmethod
    def _create_optimization(self):
        pass

    def _round(self, x, tol=1e-5):
        if abs(x - 1) < tol:
            return 1
        if abs(x) < tol:
            return 0

        return x

    def n_vars(self):
        return self.cnf_handler.n_vars

    def m_clauses(self):
        return self.cnf_handler.m_clauses

    def clauses(self):
        return self.cnf_handler.clauses

    def round(self, x):
        return [self._round(xi) for xi in x]

    def solve(self):
        s = self.solver.Solve()

        if s == 2:
            print("UNFEASIBLE")
            self.solver.Clear()
            return s, [], []

        result = self.solver.Objective().Value()
        witness = self.round([v.solution_value() if not isinstance(v, int) else v for v in self.vars])

        return s, result, witness

    def verify(self, witness):
        clauses = self.clauses()
        witness = np.array([-1 if xi == 0 else 1 for xi in witness])

        for c in clauses:
            abs_c = np.abs(c) 
            idx = abs_c - 1
            sgn = c/abs_c
            res = np.max(witness[idx]*sgn)
            if res != 1:
                return False

        return True

    def verify_partial(self, partial_witness):

        satisfied_clauses = []
        clauses = self.clauses()

        witness = np.array(partial_witness)

        for i in range(len(partial_witness)):

            if not isinstance(partial_witness[i], int):
                witness[i] = 0

            elif partial_witness[i] == 0:
                witness[i] = -1

        for i, c in enumerate(clauses):

            abs_c = np.abs(c) 
            idx = abs_c - 1
            sgn = c/abs_c
            res = np.max(witness[idx]*sgn)

            if res == 1:
                satisfied_clauses.append(i)

        return satisfied_clauses

    def load_cnf(self, filename):
        self.cnf_handler.load(filename)

    def create_lp(self):
        self._init_objects()
        self._create_optimization()
