from satlp.linear_solver.baseclass_implementation import SATasLP
from scipy.optimize import linprog
import numpy as np

# for simplex use method='highs-ds' for IPM use method='highs-ipm'

class SATasLPFeasibility(SATasLP):

    def __init__(self, filename=None, cnf_handler=None, fixing={}, method='highs-ds'):
        super().__init__(filename, cnf_handler, method)
        self.fixing = fixing

    def _init_objects(self):
        
        n_vars = self.cnf_handler.n_vars
        m_clauses = self.cnf_handler.m_clauses
        clauses = self.cnf_handler.clauses

        y_lb = np.zeros(shape=m_clauses)
        A_lb = np.zeros(shape=(m_clauses, n_vars))
        for i, c in enumerate(clauses):
            
            res_fixed = np.array([self.g(xi) for xi in c if abs(xi) in self.fixing.keys()])
            res = np.array([np.sign(xi) for xi in c if abs(xi) not in self.fixing.keys()])

            y_lb[i] = 1 - sum(res_fixed) - sum(res < 0)
            for j in c:
                idx = abs(j)-1
                if abs(j) not in self.fixing.keys():
                    A_lb[i][idx] = np.sign(j).item()

        # feasibility
        c = np.zeros(n_vars)

        # scipy linprog deals with only minimization of upperbounded matrices 
        self.y_ub = -y_lb
        self.A_ub = -A_lb
        self.c = c

    def _create_optimization(self):
        self.bounds = [0,1]

    def switch(self):
        pass

    def set_coefs_for_projection(self, decision, positive=True):
        n_vars = self.cnf_handler.n_vars
        c = np.zeros(n_vars)
        c[decision] = 1 if positive else -1
        self.c = c


class SATasLPOptimization(SATasLP):

    def __init__(
        self, 
        filename=None, 
        cnf_handler=None, 
        fixing={}, 
        last_coefs=None, 
        last_witness=None,
        method='highs-ipm',
    ):
        super().__init__(filename, cnf_handler, method)
        self.fixing = fixing
        self.last_coefs = last_coefs
        self.n_vars = cnf_handler.n_vars
        self.m_clauses = cnf_handler.m_clauses
        self.last_witness = last_witness

    def _init_objects(self):
        
        n_vars = self.cnf_handler.n_vars
        m_clauses = self.cnf_handler.m_clauses
        clauses = self.cnf_handler.clauses

        # Ax <= y
        y_lb = np.zeros(shape=m_clauses)
        A_lb = np.zeros(shape=(m_clauses, 3*n_vars))

        y_eq = np.zeros(shape=2*n_vars)
        A_eq = np.zeros(shape=(2*n_vars, 3*n_vars))

        # construct matrices following Ax >= y
        for i, c in enumerate(clauses):

            res_fixed = np.array([self.g(xi) for xi in c if abs(xi) in self.fixing.keys()])
            res = np.array([np.sign(xi) for xi in c if abs(xi) not in self.fixing.keys()])
            y_lb[i] = 1 - sum(res_fixed) - sum(res < 0)

            for j in c:
                idx = abs(j)-1
                if abs(j) not in self.fixing.keys():
                    A_lb[i][idx] = np.sign(j).item()


        # flip signs to get Ax <= y
        y_ub = -y_lb
        A_ub = -A_lb

        # optimization
        for i in range(n_vars):
            if i+1 not in self.fixing.keys():
                # x_i - y_i+ + y_i- = 1/2
                A_eq[i][i] = 1
                A_eq[i][n_vars+i] = -1
                A_eq[i][2*n_vars+i] = 1
                y_eq[i] = 1/2

            # y_n+i + y_2n+i = 1/2
                A_eq[n_vars+i][n_vars+i] = 1
                A_eq[n_vars+i][2*n_vars+i] = 1
                y_eq[n_vars+i] = 1/2

        # min x+x- + ...
        c = np.zeros(3*n_vars)
        if self.last_witness is not None:

            # get boolean array of which variable is boolean
            is_boolean = self.is_boolean(self.last_witness)
            for i in range(n_vars):

                if not is_boolean[i] and i+1 not in self.fixing.keys():
                    if self.last_witness[n_vars+i] < self.last_witness[2*n_vars+i]:
                        # x+ < x-
                        c[2*n_vars+i] = 1
                    elif self.last_witness[n_vars+i] > self.last_witness[2*n_vars+i]:
                        # x- < x+
                        c[n_vars+i] = 1
                    elif self.last_witness[n_vars+i] == self.last_witness[2*n_vars+i]:
                        c[n_vars+i] = self.last_coefs[2*n_vars+i]
                        c[2*n_vars+i] = self.last_coefs[n_vars+i]
                    # need to flip the coef if they are equal, but need last coef for this
                        
                elif is_boolean[i] and i+1 not in self.fixing.keys():
                    if self.is_one(self.last_witness[i]):
                        c[2*n_vars+i] = 1

                    elif self.is_zero(self.last_witness[i]):
                        c[n_vars+i] = 1
                    
        else:
            c[2*n_vars:] = 1

        # scipy linprog deals with only minimization of upperbounded matrices 
        self.y_ub = y_ub
        self.A_ub = A_ub
        self.y_eq = y_eq
        self.A_eq = A_eq
        self.c = c

    def _create_optimization(self):

        n = self.n_vars
        self.bounds = [
            [0,1] if i < n 
            else [0, 1/2] 
            for i in range(3*n)
        ]

    def switch(self):
        n_vars = self.cnf_handler.n_vars
        for i in range(n_vars):
            if self.c[n_vars+i] == 1/2:
                self.c[n_vars+i] = 0
                self.c[2*n_vars+i] = 1/2

            elif self.c[2*n_vars+i] == 1/2:
                self.c[2*n_vars+i] = 0
                self.c[n_vars+i] = 1/2

    def coefs_to_witness(self, coefs):
        
        n_vars = self.cnf_handler.n_vars
        witness = []
        for i in range(n_vars):
            if coefs[n_vars+i] == 1/2:
                witness.append(i+1)
            elif coefs[2*n_vars+i] == 1/2:
                witness.append(-i-1)

        return np.array(witness)

    def set_coefs_for_projection(self, decision, positive=True):

        c = np.zeros(3*n_vars)
        if positive:
            c[n_vars+decision] = 1/2

        else:
            c[2*n_vars+decision] = 1/2
        self.c = c

class SATasLPOptimizationDual(SATasLP):

    def __init__(self, filename=None, cnf_handler=None, fixing={}, method='highs-ipm'):
        super().__init__(filename, cnf_handler, method)
        self.fixing = fixing
        self.n_vars = cnf_handler.n_vars
        self.m_clauses = cnf_handler.m_clauses
        self.c = None

    def _init_objects(self):
        
        n_vars = self.cnf_handler.n_vars
        m_clauses = self.cnf_handler.m_clauses
        clauses = self.cnf_handler.clauses

        # Ax <= y
        y_lb = np.zeros(shape=m_clauses)
        A_lb = np.zeros(shape=(m_clauses, n_vars))

        # construct matrices following Ax >= y
        for i, c in enumerate(clauses):
            
            res_fixed = np.array([self.g(xi) for xi in c if abs(xi) in self.fixing.keys()])
            res = np.array([np.sign(xi) for xi in c if abs(xi) not in self.fixing.keys()])
            y_lb[i] = 1 - sum(res_fixed) - sum(res < 0)

            for j in c:
                idx = abs(j)-1
                if abs(j) not in self.fixing.keys():
                    A_lb[i][idx] = np.sign(j).item()

        # flip signs to get Ax <= y
        y_ub = -y_lb
        A_ub = -A_lb

        # optimization
        c = np.array([1.0 if i+1 not in self.fixing.keys() else 0 for i in range(n_vars)])
        self.c = c
        
        # scipy linprog deals with only minimization of upperbounded matrices 
        self.y_ub = y_ub
        self.A_ub = A_ub

    def _create_optimization(self):

        n = self.n_vars
        self.bounds = [
            [0,1]
            for i in range(n)
        ]

    def switch(self):
        self.c = -self.c

    def set_coefs_for_projection(self, idx, direction=0):
        n_vars = self.cnf_handler.n_vars
        c = np.zeros(n_vars)
        c[idx] = 1 if direction == 0 else -1
        self.c = c
