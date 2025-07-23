from numpy import array as np_array
from numpy import append as np_append

class CNFLoader():

    def __init__(self, filename=None):
        self.n_vars = None
        self.m_clauses = None
        self.clauses = None
        self.learnt_clauses = 0
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
        self.clauses = clauses

    def load(self, filename):
        if (self.clauses == None) or (filename != self.filename):
            try:
                self._load_from_file(filename)
            except FileError:
                print("Something went wrong loading the file. Make sure it is in DIMCAS formats.")

        return inf_assignments        

    def add_clause(self, clause):
        self.clauses.append(clause)
        self.m_clauses += 1