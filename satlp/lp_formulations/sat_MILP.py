from satlp.baseclass_implementation import SATasMILP

import numpy as np

class SATasMILPOptimization(SATasMILP):

    def __init__(self, filename=None, relaxed_vars=[], M=1):
        super().__init__(filename)
        self.relaxed_vars = relaxed_vars
        self.M = M

    def _init_objects(self):

        self.vars = [self.solver.NumVar(0, 1.0, f"x_{i}") for i in range(1, self.n_vars() + 1)]
        self.vars_prime = [self.solver.NumVar(0, 0.5, f"y_{i}") for i in range(1, self.n_vars() + 1)]
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

class SATasMILPFeasibility(SATasMILP):

    def __init__(self, filename=None, relaxed_vars=[]):
        super().__init__(filename)
        self.relaxed_vars = relaxed_vars

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

