from ortools.linear_solver import pywraplp
from cnf_loader import CNFLoader
import numpy as np
from sat_as_lp import SATasLP

class SATasLPAbs(SATasLP):

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


if __name__ == "__main__":

    relaxed_vars = [1,2,3,4,5,6,7,8,9,10,11,12,15,16,17,19,20]
    lp_obj = SATasLPAbs(M=1, relaxed_vars=relaxed_vars)
    # filename = "cnfs/cnf_sat.txt"
    filename = "cnfs/uf20-018.cnf"
    lp_obj.create_lp(filename)

    res, witness = lp_obj.solve()
    betas = [v.solution_value() for v in lp_obj.betas]

    print(betas)

    witness_str = {lp_obj.vars[i].name(): witness[i] for i in range(len(lp_obj.vars))}
    print(witness_str)
    print("----------------")
    
    result = lp_obj.solver.Objective().Value()
    print(result)
    print("----------------")