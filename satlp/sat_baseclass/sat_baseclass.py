from satlp import CNFLoader

from ortools.linear_solver import pywraplp
from scipy.optimize import linprog
import numpy as np
from abc import ABC, abstractmethod

class SATasLPBaseclass(ABC):

    def __init__(self, filename=None):
        self.cnf_handler = CNFLoader(filename)

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

    def add_clause(self, clause):
        self.cnf_handler.clauses.append(clause)
        self.cnf_handler.m_clauses += 1

    def solve(self):
        pass

    def verify(self, witness):
        clauses = self.clauses()
        witness = np.array([-1 if xi == 0 else 1 if xi == 1 else 0 for xi in witness])

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

    def create_lp(self, filename=None):
        if filename:
            self.load_cnf(filename)
        self._init_objects()
        self._create_optimization()