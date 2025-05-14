import numpy as np
class CNFLoader():

    def __init__(self, filename=None):
        self.n_vars = None
        self.m_clauses = None
        self.clauses = None
        if filename:
            self._load_from_file(filename)

    def _load_from_file(self, filename):

        with open(filename) as f:
            # eliminate comments, last lines and empty lines
            ignore_line = ['0', '%', 'c', '\n']
            cnf_data = [line.rstrip() for line in f if line[0] not in ignore_line]
            m_clauses = int(cnf_data[0].split()[-1])
            n_variables = int(cnf_data[0].split()[-2])
            clauses = [
                sorted([int(i) for i in l.split()[:-1]], key=abs)
                    for l in cnf_data[1:] if len(l) > 0
            ]

        self.n_vars = n_variables
        self.m_clauses = len(clauses)
        # self.clauses = np_array(clauses)
        self.clauses = clauses

    def load(self, filename):
        if (self.clauses == None) or (filename != self.filename):
            try:
                self._load_from_file(filename)
            except FileError:
                print("Something went wrong loading the file. Make sure it is in DIMCAS formats.")

    def simplify(self, assignments):
        simple_clauses = []
        for c in self.clauses:
            x = np.in1d(assignments, c)
            if not x.any():
                not_x = np.in1d(-assignments, c)
                if not_x.any():
                    relevant = -assignments[not_x]
                    simple_c = [c[i] for i in range(len(c)) if c[i] not in relevant]
                    simple_clauses.append(simple_c)                    

        return simple_clauses

    def simplify_pair(self, p, s_clauses):
        inf_assignments = []
        for c in s_clauses:
            if len(c) == 2:
                if c != p and np.in1d(np.abs(c), np.abs(p)).sum() == 2:
                    _c = np.array(c)
                    _p = np.array(p)
                    idx = (_c == _p)
                    inf_assignments.append(_p[idx].item())

        return inf_assignments        