from ortools.linear_solver import pywraplp
from cnf_loader import CNFLoader
import numpy as np
from sat_as_lp import SATasLP

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

if __name__ == "__main__":

    fixing = {5: 1, 9: 0, 1: 0}
    # filename = "cnfs/cnf_sat.txt"
    filename = "cnfs/uf20-0999.cnf"
    lp_obj = SATasLPWithFixing(filename=filename, fixing=fixing)
    # filename = "cnfs/uf20-0999.cnf"
    lp_obj.create_lp()
    status, res, witness = lp_obj.solve()

    print(lp_obj.verify(witness))
    print("----------------")
    witness_str = {lp_obj.vars[i].name(): witness[i] for i in range(len(lp_obj.vars))}
    print(witness_str)
    print("----------------")
    result = lp_obj.solver.Objective().Value()
    print(result)
    print("----------------")
    sat_clauses = lp_obj.verify_partial(witness)
    print(sat_clauses, len(sat_clauses))