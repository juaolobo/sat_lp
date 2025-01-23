from sat_baseclass.sat_as_lp import SATasLP

from ortools.linear_solver import pywraplp
import numpy as np

class SATasLPOriginal(SATasLP):

    def __init__(self, filename=None, fixing={}):
        self.fixing = fixing
        super().__init__(filename, None)

    def _init_objects(self):

        self.vars = [
            self.solver.NumVar(0 ,1.0, f"x_{i}") if i in self.relaxed_vars
            else 
                self.solver.NumVar(self.fixing[i], self.fixing[i], f"X_{i}")
            for i in range(1, self.n_vars() + 1)        
        ]
        
        self.pos_vars = [self.solver.NumVar(0.0, 0.5, f"x+_{i}") for i in range(1,self.n_vars()+1)]
        self.neg_vars = [self.solver.NumVar(0.0, 0.5, f"x-_{i}") for i in range(1,self.n_vars()+1)]

        # adjust constant terms of the inequalities
        clauses = self.clauses()
        res = [1 for c in clauses]
        for i in range(self.m_clauses()):
            res[i] = np.sum([0 if x > 0 else -1 for x in clauses[i]]) + 1

        # create clause constraints
        for i, c in enumerate(clauses):
            coefs = []
            for l in c:
                l_v = np.abs(l)
                sgn = l_v/l
                idx = l_v - 1
                coefs.append(sgn*self.vars[idx])
            
            self.solver.Add(self.solver.Sum(coefs) >= res[i])

        for v in range(self.n_vars()) :
            self.solver.Add(self.vars[v] + self.pos_vars[v] - self.neg_vars[v] == 0.5)
            self.solver.Add(self.pos_vars[v] + self.neg_vars[v] <= 0.5)
            self.solver.Add(self.pos_vars[v] + self.neg_vars[v] >= 0)        

    def _create_optimization(self):
        # no optimization function
        self.solver.Maximize(
            self.solver.Sum(
                [
                    self.solver.Sum(self.pos_vars),
                    self.solver.Sum(self.neg_vars),
                ]
            )
        )

    def create_lp(self, filename=None):

        if filename:
            self.load_cnf(filename)

        self.relaxed_vars = [
            i for i in range(1, self.n_vars()+1) 
                if i not in self.fixing.keys()
            ]
        self._init_objects()
        self._create_optimization()


class SATasLPWithFixing(SATasLP):

    def __init__(self, filename=None, fixing={}):
        self.fixing = fixing
        super().__init__(filename, None)

    def _init_objects(self):

        self.vars = [
            self.solver.NumVar(0 ,1.0, f"x_{i}") if i in self.relaxed_vars
            else 
                self.solver.NumVar(self.fixing[i], self.fixing[i], f"B_{i}")
            for i in range(1, self.n_vars() + 1)        
        ]
        
        # adjust constant terms of the inequalities
        clauses = self.clauses()
        res = [1 for c in clauses]
        for i in range(self.m_clauses()):
            res[i] = np.sum([0 if x > 0 else -1 for x in clauses[i]]) + 1

        # create clause constraints
        for i, c in enumerate(clauses):
            coefs = []
            for l in c:
                l_v = np.abs(l)
                sgn = l_v/l
                idx = l_v - 1
                coefs.append(sgn*self.vars[idx])
            
            self.solver.Add(self.solver.Sum(coefs) >= res[i])

    def _create_optimization(self):
        # no optimization function
        self.solver.Maximize(1)

    def create_lp(self, filename=None):

        if filename:
            self.load_cnf(filename)

        self.relaxed_vars = [
            i for i in range(1, self.n_vars()+1) 
                if i not in self.fixing.keys()
            ]
        self._init_objects()
        self._create_optimization()

class SATasMILP(SATasLP):

    def __init__(self, filename=None, relaxed_vars=[], M=1):
        super().__init__(filename, relaxed_vars)
        self.M = M

    def _init_objects(self):

        self.vars = [self.solver.NumVar(0, 1.0, f"x_{i}") for i in range(1, self.n_vars() + 1)]
        self.vars_prime = [self.solver.NumVar(0, 0.5, f"x_{i}") for i in range(1, self.n_vars() + 1)]
        self.betas = [
            self.solver.BoolVar(f"b_{i}") if i not in self.relaxed_vars
            else self.solver.NumVar(0 ,1.0, f"beta_{i}")
            for i in range(1, self.n_vars() + 1)        
        ]

        # adjust constant terms of the inequalities
        clauses = self.clauses()
        res = [1 for c in clauses]
        for i in range(self.m_clauses()):
            res[i] = np.sum([0 if x > 0 else -1 for x in clauses[i]]) + 1

        # create clause constraints
        for i, c in enumerate(clauses):
            coefs = []
            for l in c:
                l_v = np.abs(l)
                sgn = l_v/l
                idx = l_v - 1
                coefs.append(sgn*self.vars[idx])
            
            self.solver.Add(self.solver.Sum(coefs) >= res[i])

        # add absolute value contraints
        for i in range(self.n_vars()):

            self.solver.Add(self.vars[i]-1/2 + self.M*self.betas[i] >= self.vars_prime[i])
            self.solver.Add(-self.vars[i]+1/2 + self.M*(1 - self.betas[i]) >= self.vars_prime[i])
            self.solver.Add(self.vars[i]-1/2 <= self.vars_prime[i])
            self.solver.Add(-self.vars[i]+1/2 <= self.vars_prime[i])

    def _create_optimization(self):
        # optimizing for sum of artificial variables
        self.solver.Maximize(
            self.solver.Sum(self.vars_prime)
        )

class SATasMILPSimple(SATasLP):

    def __init__(self, filename=None, relaxed_vars=[]):
        super().__init__(filename, relaxed_vars)

    def _init_objects(self):

        self.vars = [
            self.solver.BoolVar(f"b_{i}") if i not in self.relaxed_vars
            else self.solver.NumVar(0 ,1.0, f"beta_{i}")
            for i in range(1, self.n_vars() + 1)        
        ]

        # adjust constant terms of the inequalities
        clauses = self.clauses()
        res = [1 for c in clauses]
        for i in range(self.m_clauses()):
            res[i] = np.sum([0 if x > 0 else -1 for x in clauses[i]]) + 1

        # create clause constraints
        for i, c in enumerate(clauses):
            coefs = []
            for l in c:
                l_v = np.abs(l)
                sgn = l_v/l
                idx = l_v - 1
                coefs.append(sgn*self.vars[idx])
            
            self.solver.Add(self.solver.Sum(coefs) >= res[i])

    def _create_optimization(self):
        # no optimization func
        self.solver.Maximize(1)

