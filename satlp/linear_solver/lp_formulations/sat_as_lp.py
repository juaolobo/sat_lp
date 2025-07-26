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


class SATasLPOptimization(SATasLP):

    def __init__(
        self, 
        filename=None, 
        cnf_handler=None, 
        fixing={}, 
        potential_coefs=None, 
        method='highs-ipm',
    ):
        super().__init__(filename, cnf_handler, method)
        self.fixing = fixing
        self.potential_coefs = potential_coefs
        self.n_vars = cnf_handler.n_vars
        self.m_clauses = cnf_handler.m_clauses

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
                # x_i - y_n+i + y_2n+i = 1/2
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
        if self.potential_coefs is not None:

            # get boolean array of which variable is boolean
            is_boolean = self.is_boolean(self.potential_coefs)
            for i in range(n_vars):

                if not is_boolean[i] and i+1 not in self.fixing.keys():
                    if self.potential_coefs[n_vars+i] < self.potential_coefs[2*n_vars+i]:
                        # c[2*n_vars+i] = self.potential_coefs[n_vars+i]
                        c[2*n_vars+i] = 1/2
                    else:
                        # c[n_vars+i] = self.potential_coefs[2*n_vars+i]
                        c[n_vars+i] = 1/2
                        
                elif is_boolean[i] and i+1 not in self.fixing.keys():
                    if self.is_one(self.potential_coefs[i]):
                        c[2*n_vars+i] = 1/2

                    elif self.is_zero(self.potential_coefs[i]):
                        c[n_vars+i] = 1/2
        else:
            c[2*n_vars:] = 1/2

    
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


class SATasLPOptimization2(SATasLP):

    def __init__(self, filename=None, cnf_handler=None, fixing={}, potential_coefs=None, method='highs-ipm'):
        super().__init__(filename, cnf_handler, method)
        self.fixing = fixing
        self.potential_coefs = potential_coefs
        self.n_vars = cnf_handler.n_vars
        self.m_clauses = cnf_handler.m_clauses

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

        y_eq = np.zeros(shape=n_vars)
        A_eq = np.zeros(shape=(n_vars, n_vars))

        # optimization
        # min x.x + ...
        if self.potential_coefs is not None:
            c = 1-2*self.potential_coefs

        else:
            c = np.zeros(n_vars)

        print(f"C: {c}")
        print(f"POT: {self.potential_coefs}")

        # scipy linprog deals with only minimization of upperbounded matrices 
        self.y_ub = y_ub
        self.A_ub = A_ub
        self.y_eq = y_eq
        self.A_eq = A_eq
        self.c = c

    def _create_optimization(self):

        n = self.n_vars
        self.bounds = [
            [0,1]
            for i in range(n)
        ]
