from satlp.linear_solver import (
    SATasLPOptimization,
    SATasLPFeasibility,
)
from satlp.boolean_solver import BooleanSolver
from line_profiler import profile

class HybridSolver:

    def __init__(self, filename, lp_solver, solutions='', track_history=False):

        self.fixing = {}
        self.lp_solver = lp_solver(fixing=self.fixing, filename=filename, method='highs-ipm')
        self.bool_solver = BooleanSolver(filename, verbose=0)
        self.track_history = track_history
        self.solution_history = []

    def check_solutions(self, solutions, solution):
        errors = []
        min_err = n_vars
        argmin = -1
        for l in solutions:
            real_solution = np.array([int(x) for x in l.split()[:-1]])
            error = sum(real_solution != solution) - sum(solution == 0)
            errors.append(error)
            if error < min_err:
                min_err = error
                min_sol = real_solution
                wrong_idxs = real_solution[(real_solution != solution) & (solution != 0)]

        return min_err, min_sol, wrong_idxs

    @profile
    def solve_linear(self):

        self.lp_solver.create_lp()
        return self.lp_solver.solve()

    @profile
    def solve_boolean(self):

        # assert fixing = self.fixing
        linear_sol = [xi if self.fixing[xi] else -xi for xi in self.fixing.keys()]

        self.bool_solver.restart()
        resolved, formula = self.bool_solver.propagate_linear(linear_sol)
        new_clauses = [f.clause for f in formula.formula[self.lp_solver.m_clauses():]]

        if len(new_clauses) == 0:
            resolved, formula = self.bool_solver.extend_solution()
            new_clauses = [f.clause for f in formula.formula[self.lp_solver.m_clauses():]]

        for c in new_clauses:
            self.lp_solver.add_clause(c)
            self.bool_solver.list_clause.append(c)

        return resolved

    @profile
    def solve(self):
        it = 0
        witness = self.solve_linear()
        while not self.lp_solver.verify(witness):
            # i.e. INFEASIBLE
            if witness == None:
                witness = [0.5 for _ in range(lp.n_vars())]
                return None

            fixing = {i+1: xi for i, xi in enumerate(witness) if xi.is_integer()}

            # track solution history
            if self.track_history:
                self.solution_history.append(fixing)

            # if we hit a fixed point of the linear solver
            if fixing == self.fixing:
                resolved = self.solve_boolean()
                self.fixing = {abs(xi): 1.0 if xi > 0 else 0.0 for xi in resolved}

            # else continue to evolve the linear solution
            else:
                self.fixing = fixing

            if it%100 == 0:
                dimacs = [xi if ki else -xi for xi, ki in self.fixing.items()]
                print(f"Current solution (iteration {it}): {dimacs}")

            self.lp_solver.restart(fixing=self.fixing)                
            witness = self.solve_linear()

            it += 1        

        print(f"Finished in {it} iterations")
        return witness

    def verify(self, witness):
        if witness:
            if self.lp_solver.verify(witness):
                print("SATISFIABLE")
                linear_sol = [i+1 if xi == 1 else -i-1 for i,xi in enumerate(witness)]
                print(f"WITNESS: {linear_sol}")
        
        else: print("UNSATISFIABLE")


