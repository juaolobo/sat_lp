from ortools.linear_solver import pywraplp
from cnf_loader import CNFLoader
import numpy as np
from sat_as_lp import SATasLP

class SATasLPSimple(SATasLP):

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


if __name__ == "__main__":

    relaxed_vars = [1,2,3,4,5,6,7,8,9,10,11,12,15,16,17,19,20]
    lp_obj = SATasLPSimple(relaxed_vars=relaxed_vars)
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