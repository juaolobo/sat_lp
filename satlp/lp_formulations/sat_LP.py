from satlp.baseclass_implementation import SATasLP

from scipy.optimize import linprog
import numpy as np

# for simplex use method='highs-dm' for IPM use method='highs-ipm'

class SATasLPFeasibility(SATasLP):

    def __init__(self, filename=None, fixing={}, method='highs-dm'):
        self.fixing = fixing
        super().__init__(filename, method)

    def _init_objects(self):
        
        y_lb = np.zeros(shape=self.m_clauses())
        A_lb = np.zeros(shape=(self.m_clauses(), self.n_vars()))
        for i, c in enumerate(self.clauses()):

            res_fixed = np.array([self.g(xi) for xi in c if abs(xi) in self.fixing.keys()])
            res = np.array([np.sign(xi) for xi in c if abs(xi) not in self.fixing.keys()])

            y_lb[i] = 1 - sum(res_fixed) - sum(res < 0)
            for j in c:
                idx = abs(j)-1
                if abs(j) not in self.fixing.keys():
                    A_lb[i][idx] = np.sign(j).item()

        # feasibility
        c = np.zeros(self.n_vars())

        # scipy linprog deals with only minimization of upperbounded matrices 
        self.y_ub = -y_lb
        self.A_ub = -A_lb
        self.c = c

    def _create_optimization(self):
        self.bounds = [0,1]

class SATasLPOptimization(SATasLP):

    def __init__(self, filename=None, fixing={}, method='highs-dm'):
        self.fixing = fixing
        super().__init__(filename, method)

    def _init_objects(self):
        
        # Ax <= y
        y_lb = np.zeros(shape=self.m_clauses())
        A_lb = np.zeros(shape=(self.m_clauses(), 3*self.n_vars()))

        y_eq = np.zeros(shape=self.n_vars())
        A_eq = np.zeros(shape=(self.n_vars(), 3*self.n_vars()))

        # construct matrices following Ax >= y
        for i, c in enumerate(self.clauses()):
            
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

        n = self.n_vars()
        y_ub1 = np.zeros(shape=self.n_vars())
        A_ub1 = np.zeros(shape=(self.n_vars(), 3*self.n_vars()))

        # y+ + y- <= 1/2
        for i in range(n):
            if i+1 not in self.fixing.keys():
                A_ub1[i][n+i] = 1
                A_ub1[i][2*n+i] = 1
                y_ub1[i] = 1/2

        # # inefficient but more readable
        y_ub = np.concatenate([y_ub, y_ub1])
        A_ub = np.concatenate([A_ub, A_ub1])
        # optimization
        # x_i + y_n+i - y_2n+i = 1/2
        for i in range(self.n_vars()):
            if i+1 not in self.fixing.keys():
                A_eq[i][i] = 1
                A_eq[i][n+i] = 1
                A_eq[i][2*n+i] = -1
                y_eq[i] = 1/2

        c = -np.ones(3*self.n_vars())
        c[:n] = 0
        for i in range(n):
            if i+1 in self.fixing.keys():
                c[i+n] = 0
                c[i+2*n] = 0

        # scipy linprog deals with only minimization of upperbounded matrices 
        self.y_ub = y_ub
        self.A_ub = A_ub
        self.y_eq = y_eq
        self.A_eq = A_eq
        self.c = c

    def _create_optimization(self):

        n = self.n_vars()
        self.bounds = [
            [0,1] if i < n 
            else [0, 1/2] 
            for i in range(3*n)
        ]