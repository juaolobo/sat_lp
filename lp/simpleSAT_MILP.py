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

    def verify(self, witness):
        clauses = self.clauses()
        witness = np.array([-1 if xi == 0 else 1 for xi in witness])
        for c in clauses:
            abs_c = np.abs(c) 
            idx = abs_c - 1
            sgn = c/abs_c
            res = np.max(witness[idx]*sgn)
            if res != 1:
                return False

        return True

    def create_lp(self, filename=None):
        self.cnf_handler.load(filename)
        self._init_objects()
        self._create_optimization()


if __name__ == "__main__":

    relaxed_vars = (1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17, 19, 20)
    # relaxed_vars = (1, 2, 3, 4, 6, 7, 8, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19)
    lp_obj = SATasLP(relaxed_vars=relaxed_vars)
    # filename = "cnfs/cnf_sat.txt"
    filename = "cnfs/uf20-018.cnf"
    # filename = "cnfs/uf20-0999.cnf"
    lp_obj.create_lp(filename)
    res, witness = lp_obj.solve()
    print(lp_obj.verify(witness))
    print("----------------")
    witness_str = {lp_obj.vars[i].name(): witness[i] for i in range(len(lp_obj.vars))}
    print(witness_str)
    print("----------------")
    result = lp_obj.solver.Objective().Value()
    print(result)
    print("----------------")