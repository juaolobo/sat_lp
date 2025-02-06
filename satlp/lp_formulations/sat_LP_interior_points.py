from satlp.baseclass_implementation import SATasLPIP

from scipy.optimize import linprog
import numpy as np

class SATasLPFeasibilityIP(SATasLPIP):

    def __init__(self, filename=None, fixing={}):
        self.fixing = fixing
        self.g = lambda x: 1 - self.fixing[abs(x)] if x < 0 else self.fixing[abs(x)]
        super().__init__(filename)

    def _init_objects(self):
        
        y = np.zeros(shape=self.m_clauses())
        A = np.zeros(shape=(self.m_clauses(), self.n_vars()))
        for i, c in enumerate(self.clauses()):
            
            res_fixed = np.array([self.g(xi) for xi in c if abs(xi) in self.fixing.keys()])
            res = np.array([np.sign(xi) for xi in c if abs(xi) not in self.fixing.keys()])

            y[i] = 1 - sum(res_fixed) - sum(res < 0)
            for j in c:
                idx = abs(j)-1
                if abs(j) not in self.fixing.keys():
                    A[i][idx] = np.sign(j).item()

        # feasibility

        c = np.zeros(self.n_vars())


        # scipy linprog deals with only minimization of upperbounded matrices 
        self.y = -y
        self.A = -A
        self.c = c

    def _create_optimization(self):
        self.bounds = [0,1]
        self.method = 'highs-ipm'
