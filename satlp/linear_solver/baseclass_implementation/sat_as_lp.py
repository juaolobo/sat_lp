from satlp.linear_solver.sat_baseclass.sat_baseclass import SATasLPBaseclass

from scipy.optimize import linprog
import numpy as np
from abc import ABC, abstractmethod

class SATasLP(SATasLPBaseclass):

    def __init__(self, filename=None, cnf_handler=None, method=None, eps=10e-6):
        super().__init__(filename, cnf_handler)
        self.eps = eps
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
        res = result.fun
        print(f"LAST OPTIMIZATION RESULT: {res}")

        x = result.x

        if result.success:
            witness = np.array(
                [
                    x[i-1].item() 
                    if i not in self.fixing.keys() 
                    else self.fixing[i] 
                    for i in range(1,len(x)+1)
                ]
            )

            return witness

        print("INFEASIBLE")


    def verify(self, witness):

        clauses = self.cnf_handler.clauses
        ver_witness = np.array([-1 if xi == 0 else 1 if xi == 1 else 0 for xi in witness])
        for c in clauses:
            abs_c = np.abs(c) 
            idx = abs_c - 1
            sgn = c/abs_c
            res = np.max(ver_witness[idx]*sgn)
            if res != 1:
                return False

        return True

    def check_conflict(self, witness):

        clauses = self.cnf_handler.clauses
        n_vars = self.cnf_handler.n_vars
        ver_witness = np.array([-1 if xi == 0 else 1 if xi == 1 else 0 for xi in witness[:n_vars]])

        for c in clauses:
            abs_c = np.abs(c) 
            idx = abs_c - 1
            sgn = c/abs_c
            res = np.max(ver_witness[idx]*sgn)

            if res == -1:
                return True

        return False

    def check_blocked(self, witness):

        clauses = self.cnf_handler.clauses
        n_vars = self.cnf_handler.n_vars
        ver_witness = np.array([-1 if xi == 0 else 1 if xi == 1 else 0 for xi in witness[:n_vars]])

        if self.verify(witness):
            return False

        for i, xi in enumerate(ver_witness):
            if xi == 0:
                # pick xi = [-1, 1] and check for conflicts
                aux = ver_witness[i]
                for v in [-1, 1]:
                    ver_witness[i] = v
                    if not self.check_conflict(ver_witness):
                        return False

                ver_witness[i] = aux

        return True

    def linear_to_witness(self, witness):
        n_vars = self.cnf_handler.n_vars
        solution = []

        for i in range(n_vars):
            if self.is_one(witness[i]):
                solution.append(i+1)
            elif self.is_zero(witness[i]):
                solution.append(-i-1)

        return solution

    def get_active_clauses(self, partial_witness):

        sat_clauses = []
        unsat_clauses = []
        clauses = self.cnf_handler.clauses

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

    def restart(self, fixing={}, last_coefs=None, last_witness=None):
        self.fixing = fixing
        self.last_coefs = last_coefs
        self.last_witness = last_witness

    def is_one(self, x):
        return np.abs(1-x) < self.eps

    def is_zero(self, x):
        return np.abs(x) < self.eps

    def is_boolean(self, x):
        return self.is_one(x) | self.is_zero(x)
