from ortools.linear_solver import pywraplp
from cnf_loader import CNFLoader
import numpy as np

class SATasLP:

    def __init__(self, filename=None, relaxed_vars=[], M=1):
        self.cnf_handler = CNFLoader(filename)
        self.relaxed_vars = relaxed_vars
        self.M = M
        self.solver = pywraplp.Solver.CreateSolver("SCIP")
        if not self.solver:
            raise Exception("Solver creation failed")

    def n_vars(self):
        return self.cnf_handler.n_vars

    def m_clauses(self):
        return self.cnf_handler.m_clauses

    def clauses(self):
        return self.cnf_handler.clauses

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

    def _round(self, x, tol=1e-5):
        if abs(x - 1) < tol:
            return 1
        if abs(x) < tol:
            return 0

        return x

    def round(self, x):
        return [self._round(xi) for xi in x]

    def solve(self):
        self.solver.Solve()
        result = self.solver.Objective().Value()
        witness = self.round([v.solution_value() for v in self.vars])
        
        return result, witness

    def create_lp(self, filename=None):
        self.cnf_handler.load(filename)
        self._init_objects()
        self._create_optimization()


if __name__ == "__main__":

    lp_obj = SATasLP()
    # filename = "cnfs/cnf_sat.txt"
    filename = "cnfs/uf20-018.cnf"
    lp_obj.create_lp(filename)
    res, witness = lp_obj.solve()
    witness_str = {lp_obj.vars[i].name(): witness[i] for i in range(len(lp_obj.vars))}
    print(witness_str)
    print("----------------")
    result = lp_obj.solver.Objective().Value()
    print(result)
    print("----------------")