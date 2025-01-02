class CNFLoader():

    def __init__(self, filename=None):
        self.n_vars = None
        self.m_clauses = None
        self.clauses = None
        if file:
            self._load_from_file(filename)


    def _load_from_file(self, filename):

        with open(filename) as f:
            cnf_data = [line.rstrip() for line in f if line != '\n']
            m_clauses = int(cnf_data[0].split()[-1])
            n_variables = int(cnf_data[0].split()[-2])
            clauses = [
                [int(i) for i in l.split()[:-1]]
                    for l in cnf_data[1:]
            ]

        self.n_vars = n_variables
        self.m_clauses = m_clauses
        self.clauses = clauses

    def load(self, filename):
        try:
            self._load_from_file(filename)
        except FileError:
            print("Something went wrong loading the file. Make sure it is in DIMCAS formats.")
