from satlp.linear_solver import (
    SATasLPOptimization,
    SATasLPFeasibility,
)
from satlp.boolean_solver import BooleanSolver
from satlp.cnf_loader import CNFLoader

class HybridSolver:

    def __init__(self, filename, lp_solver, method='highs-ipm', solutions='', track_history=False):

        self.fixing = {}
        self.cnf_handler = CNFLoader(filename)
        self.lp_solver = lp_solver(fixing=self.fixing, filename=filename, method=method, cnf_handler=self.cnf_handler)
        self.bool_solver = BooleanSolver(filename, verbose=0, cnf_handler=self.cnf_handler)
        self.track_history = track_history
        self.solution_history = []
        self.linear_it = 0
        self.boolean_it = 0

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

    def solve_linear(self):

        self.lp_solver.create_lp()
        return self.lp_solver.solve()

    def solve_boolean(self):
        
        # assert fixing = self.fixing
        linear_sol = [xi if self.fixing[xi] else -xi for xi in self.fixing.keys()]

        resolved, formula = self.bool_solver.propagate_linear(linear_sol)

        # UNSAT
        if formula is None:
            return None

        new_clauses = [f.clause for f in formula.formula[self.lp_solver.m_clauses():]]

        if len(new_clauses) == 0:
            resolved, formula = self.bool_solver.extend_solution()

            # UNSAT
            if formula is None:
                return None

            new_clauses = [f.clause for f in formula.formula[self.lp_solver.m_clauses():]]

        for c in new_clauses:
            self.cnf_handler.add_clause(c)
            self.cnf_handler.learnt_clauses += 1

        self.bool_solver.restart()

        return resolved

    def solve(self):
        it = 0
        witness = self.solve_linear()
        self.linear_it += 1
        while not self.lp_solver.verify(witness):
            # i.e. INFEASIBLE
            if witness == None:
                witness = [0.5 for _ in range(lp.n_vars())]
                return None

            fixing = {i+1: xi for i, xi in enumerate(witness) if xi.is_integer()}

            # if we hit a fixed point of the linear solver
            if fixing == self.fixing:
                # track solution history
                if self.track_history:
                    sol = tuple([
                            xi if self.fixing[xi] 
                            else -xi 
                            for xi in list(self.fixing.keys()) 
                            if abs(xi) <= self.lp_solver.n_vars()
                        ]
                    )
                    if len(sol) > 0:
                        self.solution_history.append(sol)
                    
                resolved = self.solve_boolean()
                self.boolean_it += 1

                # UNSAT
                if resolved is None:
                    return None
                                                           
                self.fixing = {abs(xi): 1.0 if xi > 0 else 0.0 for xi in resolved}

            # else continue to evolve the linear solution
            else:
                self.fixing = fixing

            if it%100 == 0:
                items = [(k,v) for k,v in self.fixing.items() if k < self.lp_solver.n_vars()]
                dimacs = [xi if vi else -xi for xi, vi in items]
                print(f"Current solution (iteration {it}): {dimacs}")
                print(f"Current number of clauses {self.lp_solver.m_clauses()}")
                print(f"----------------------------------------------------")

                self.fixing = {}

            self.lp_solver.restart(fixing=self.fixing)                
            witness = self.solve_linear()
            self.linear_it += 1


            it += 1        

        print(f"Finished in {it} iterations")
        return witness

    def new_solver(self):
        it = 0
        witness = self.solve_linear()
        self.linear_it += 1
        while not self.lp_solver.verify(witness):
            # i.e. INFEASIBLE
            if witness == None:
                witness = [0.5 for _ in range(lp.n_vars())]
                return None

            fixing = {i+1: xi for i, xi in enumerate(witness) if xi.is_integer()}

            if fixing == self.fixing:
                # satisfied clauses create a conflict with every single not satisfied on the boolean domain
                sat_idx, unsat_idx = self.lp_solver.get_active_clauses(witness)
                # function to check conflict/inexpansionability
                conflict = self.check_linear_conflict(witness)
                # track solution history

                linear_sol = [xi if self.fixing[xi] else -xi for xi in self.fixing.keys()]
                for c in self.bool_solver.formula.formula[unsat_idx]:
                    # pick a variable to satisfy the clause
                    # resolve conflict
                    # learn clauses
                    new_clauses = self.bool_solver.expand_and_learn(linear_sol, c)
                    print(new_clauses)
                    # if no new clauses are generated, that means that the solution is still expansionable
                    for c in new_clauses:
                        self.cnf_handler.add_clause(c)
                        self.cnf_handler.learnt_clauses += 1
                    
                    self.bool_solver.restart()
                    print(len(self.bool_solver.formula.formula))

                breakpoint()
                print(self.lp_solver.m_clauses())
                # resolved = self.solve_boolean()
                self.boolean_it += 1
                # UNSAT
                resolved = []
                if resolved is None:
                    return None

                self.fixing = {}

            # else continue to evolve the linear solution
            else:
                self.fixing = fixing

            self.lp_solver.restart(fixing=self.fixing)                
            witness = self.solve_linear()
            self.linear_it += 1


    def verify(self, witness):
        
        sat = False
        if witness:
            sat = self.lp_solver.verify(witness)
            if sat:
                print("SATISFIABLE")
                witness = self.lp_solver.linear_to_witness(witness)
                print(f"WITNESS: {witness[:self.lp_solver.n_vars()]}")
                        
            else: print("UNSATISFIABLE")

        return sat

    def check_linear_conflict(self, witness):
        conflict = True
        for i in range(self.lp_solver.n_vars()):
            if not witness[i].is_integer():
                aux = witness[i]
                for w in [0,1]:
                    witness[i] = w
                    if (self.lp_solver.A_ub @ witness <= self.lp_solver.y_ub).all():
                        conflict = False
                        break
                witness[i] = aux

        return conflict

