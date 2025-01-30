from satlp.baseclass_implementation import SATasLPIP

from scipy.optimize import linprog
import numpy as np

class SATasLPFeasibilityIP(SATasLPIP):

    def __init__(self, filename=None, fixing={}):
        self.fixing = fixing
        super().__init__(filename)

    def _init_objects(self):
        
        y = np.zeros(shape=self.m_clauses())
        A = np.zeros(shape=(self.m_clauses(), self.n_vars()))

        for i, c in enumerate(self.clauses()):
            y[i] = 1 - sum(np.sign(c) < 0)
            for j in c:
                idx = abs(j)-1
                A[i][idx] = np.sign(j).item() if j not in self.fixing.keys() else self.fixing[j]

        # feasibility

        c = np.zeros(self.n_vars())

        # scipy linprog deals with only minimization of upperbounded matrices 
        self.y = -y
        self.A = -A
        self.c = c

    def _create_optimization(self):
        self.bounds = [0,1]
        self.method = 'highs'
