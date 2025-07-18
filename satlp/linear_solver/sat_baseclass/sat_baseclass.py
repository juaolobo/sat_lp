from satlp.cnf_loader import CNFLoader

from scipy.optimize import linprog
import numpy as np
from abc import ABC, abstractmethod

class SATasLPBaseclass(ABC):

    def __init__(self, filename=None, cnf_handler=None):
        self.cnf_handler = cnf_handler if cnf_handler is not None else CNFLoader(filename)

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
        pass

    def verify(self, witness, translate=False):

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

    def check_conflict(self, witness, translate=False):

        clauses = self.clauses()
        witness = np.array([-1 if xi == 0 else 1 if xi == 1 else 0 for xi in witness])
        for c in clauses:
            abs_c = np.abs(c) 
            idx = abs_c - 1
            sgn = c/abs_c
            res = (witness[idx]*sgn == -1).all()
            print(witness[idx]*sgn)
            if res:
                return True

        return False

    def linear_to_witness(self, witness):
        solution = [i+1 if xi == 1 else -i-1 for i,xi in enumerate(witness)]
        return solution

    def get_active_clauses(self, partial_witness):

        sat_clauses = []
        unsat_clauses = []
        clauses = self.clauses()

        witness = np.array([-1 if xi == 0 else 1 if xi == 1 else 0 for xi in partial_witness])

        for i, c in enumerate(clauses):

            abs_c = np.abs(c) 
            idx = abs_c - 1
            sgn = c/abs_c
            res = np.max(witness[idx]*sgn)

            if res == 1:
                sat_clauses.append(i)
            else:
                unsat_clauses.append(i)

        return sat_clauses, unsat_clauses

    def load_cnf(self, filename):
        self.cnf_handler.load(filename)

    def create_lp(self, filename=None):
        if filename:
            self.load_cnf(filename)
        self._init_objects()
        self._create_optimization()

    def restart(self, fixing={}):
        self.fixing = fixing
        self._init_objects()
        self._create_optimization()
