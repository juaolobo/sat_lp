from satlp.linear_solver import (
    SATasLPOptimization,
    SATasLPFeasibility,
)
from satlp.boolean_solver import BooleanSolver
from satlp.cnf_loader import CNFLoader
import numpy as np

class HybridSolver:

    def __init__(
        self,
        filename,
        lp_solver,
        method='highs-ipm',
        solutions='',
        track_history=False
    ):

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

        resolved, new_clauses = self.bool_solver.propagate_linear(linear_sol)

        # UNSAT
        if new_clauses is None:
            return None

        if len(new_clauses) == 0:
            resolved, new_clauses = self.bool_solver.extend_solution()

            # UNSAT
            if new_clauses is None:
                return None

        for c in new_clauses:
            self.cnf_handler.add_clause(c)
            self.cnf_handler.learnt_clauses += 1

        self.bool_solver.restart()

        return resolved

    def solve(self):

        it = 0
        n_vars = self.cnf_handler.n_vars
        m_clauses = self.cnf_handler.m_clauses
        witness = self.solve_linear()
        self.linear_it += 1

        while not self.lp_solver.verify(witness):
            # i.e. INFEASIBLE
            if witness is None:
                return None

            fixing = {
                        i+1: xi for i, xi in enumerate(witness) 
                            if xi.is_integer() 
                            if i < n_vars
                    }

            # if we hit a fixed point of the linear solver
            if fixing == self.fixing:
                # track solution history
                if self.track_history:
                    sol = tuple([
                            xi if self.fixing[xi] 
                            else -xi 
                            for xi in list(self.fixing.keys()) 
                            if xi <= n_vars
                        ]
                    )
                    if len(sol) > 0:
                        self.solution_history.append(sol)
                    
                resolved = self.solve_boolean()
                self.boolean_it += 1

                # UNSAT
                if resolved is None:
                    print("UNSATISFIABLE")
                    return None
                                                           
                self.fixing = {abs(xi): 1.0 if xi > 0 else 0.0 for xi in resolved}

            # else continue to evolve the linear solution
            else:
                self.fixing = fixing

            if it%100 == 0:
                items = [(k,v) for k,v in self.fixing.items() if k < n_vars]
                dimacs = [xi if vi else -xi for xi, vi in items]
                print(f"Current solution (iteration {it}): {dimacs}")
                print(f"Current number of clauses {m_clauses}")
                print(f"----------------------------------------------------")

                self.fixing = {}

            self.lp_solver.restart(fixing=self.fixing)
            witness = self.solve_linear()
            self.linear_it += 1

            it += 1        

        print(f"Finished in {it} iterations")
        return witness

    def optimize1(self, verbose=False):

        it = 0
        witness = self.solve_linear()
        n_vars = self.cnf_handler.n_vars
        m_clauses = self.cnf_handler.m_clauses
        stops = 0
        solution_history = set()
        tracking = 0

        if verbose:
            print(f"X: {witness[:n_vars]}")
            print(f"X_+: {witness[n_vars:2*n_vars]}")
            print(f"C_+: {self.lp_solver.c[n_vars:2*n_vars]}")
            print(f"X_-: {witness[2*n_vars:]}")
            print(f"C_-: {self.lp_solver.c[2*n_vars:]}")

        dmacs = set()
        for i, xi in enumerate(witness[:n_vars]):
            if xi == 1:
                dmacs = dmacs | {i+1}
            elif xi == 0:
                dmacs = dmacs | {-i-1}

        self.lp_solver.restart(potential_coefs=witness)
        while not self.lp_solver.verify(witness):
            
            last = dmacs
            witness = self.solve_linear()

            # INFEASBIEL
            if witness is None:
                return None

            dmacs = set()
            for i, xi in enumerate(witness[:n_vars]):
                if xi == 1:
                    dmacs.add(i+1)
                elif xi == 0:
                    dmacs.add(-i-1)

            solution_history.add(tuple(self.lp_solver.c))
            tracking += 1
            self.linear_it += 1

            if verbose:
                print("------------------------------------------------------------")
                print(f"NUMBER OF DIFFERENT SOLUTIONS: {len(solution_history)}, NUMBER OF CURRENT CYCLE ITERATIONS: {tracking}")

                print(f"X: {witness[:n_vars]}")
                print(f"X_+: {witness[n_vars:2*n_vars]}")
                print(f"C_+: {self.lp_solver.c[n_vars:2*n_vars]}")
                print(f"X_-: {witness[2*n_vars:]}")
                print(f"C_-: {self.lp_solver.c[2*n_vars:]}")

                print(f"WITNESS: {sorted(list(dmacs), key=abs)}")


            if len(solution_history) < tracking:
                sat_idx, unsat_idx = self.lp_solver.get_active_clauses(witness)
                # unsat_idx = np.random.choice(unsat_idx, 1)
                linear_sol = dmacs

                new_clauses = []
                self.bool_solver.restart()
                for idx in unsat_idx:
                    # pick a variable to satisfy the clause
                    # resolve conflict
                    # learn clauses
                    learned = self.bool_solver.expand_and_learn(linear_sol, idx)

                    # UNSAT
                    if learned is None:
                        return None

                    new_clauses += learned
                    # if no new clauses are generated, that means that the solution is still expansionable

                new_clauses = set([tuple(sorted(c, key=abs)) for c in new_clauses])

                for c in new_clauses:
                    self.cnf_handler.add_clause(c)
                    self.cnf_handler.learnt_clauses += 1
                
                solution_history = set()
                tracking = 0
                self.boolean_it += 1
                self.lp_solver.restart()

            else:
                self.lp_solver.restart(potential_coefs=witness)

            self.lp_solver.create_lp()

        print("------------------------------------------------------------")
        print(f"NUMBER OF DIFFERENT SOLUTIONS: {len(solution_history)}, NUMBER OF ITERATIONS: {self.linear_it}")


        return witness      


    def optimize2(self, verbose=False):

        it = 0
        n_vars = self.cnf_handler.n_vars
        m_clauses = self.cnf_handler.m_clauses
        solution_history = set()
        tracking = 0
        test_c = set()
        witness = self.solve_linear()
        while not self.lp_solver.verify(witness):
            
            fixing = {
                i+1: xi for i, xi in enumerate(witness) 
                    if xi.is_integer() 
                    if i < n_vars
            }

            dmacs = set()
            for i, xi in enumerate(witness[:n_vars]):
                if xi == 1:
                    dmacs.add(i+1)
                elif xi == 0:
                    dmacs.add(-i-1)

            solution_history.add(tuple(witness))
            tracking += 1
            last_witness = witness
            self.linear_it += 1

            if verbose:
                print("------------------------------------------------------------")
                print(f"NUMBER OF DIFFERENT SOLUTIONS: {len(solution_history)}, NUMBER OF CURRENT CYCLE ITERATIONS: {tracking}")

                print(f"X: {witness[:n_vars]}")
                print(f"X_+: {witness[n_vars:2*n_vars]}")
                print(f"C_+: {self.lp_solver.c[n_vars:2*n_vars]}")
                print(f"X_-: {witness[2*n_vars:]}")
                print(f"C_-: {self.lp_solver.c[2*n_vars:]}")

                print(f"WITNESS: {sorted(list(dmacs), key=abs)}")


            if len(solution_history) < tracking:
                # satisfied clauses create a conflict with every single not satisfied on the boolean domain
                sat_idx, unsat_idx = self.lp_solver.get_active_clauses(witness)
                linear_sol = [xi if self.fixing[xi] else -xi for xi in self.fixing.keys()]
                # unsat_idx = np.random.choice(unsat_idx, 1, replace=False)

                self.bool_solver.restart()
                new_clauses = []
                for idx in unsat_idx:
                    # pick a variable to satisfy the clause
                    # resolve conflict
                    # learn clauses
                    learned = self.bool_solver.expand_and_learn(linear_sol, idx)
                    # UNSAT
                    if learned is None:
                        return None
                    new_clauses += learned
                    # if no new clauses are generated, that means that the solution is still expansionable
                    self.bool_solver.restart()

                new_clauses = set([tuple(c) for c in new_clauses])
                for c in new_clauses:
                    self.cnf_handler.add_clause(c)
                    self.cnf_handler.learnt_clauses += 1


                solution_history = set()
                tracking = 0
                self.boolean_it += 1
                self.fixing = dict()
                self.lp_solver.restart(fixing=self.fixing, potential_coefs=None)

            # else continue to evolve the linear solution
            else:
                self.fixing = fixing
                self.lp_solver.restart(fixing=self.fixing, potential_coefs=witness)

            witness = self.solve_linear()

            # INFEASIBLE
            if witness is None:
                return None

            self.linear_it += 1

        return witness

    def verify(self, witness):
        
        sat = False
        n_vars = self.cnf_handler.n_vars
        if witness is not None:
            sat = self.lp_solver.verify(witness)
            if sat is True:
                print("SATISFIABLE")
                witness = self.lp_solver.linear_to_witness(witness)
                print(f"WITNESS: {witness[:n_vars]}")

            else:
                print("SOMETHING IS WRONG")

            return sat

        print("UNSATISFIABLE")
        return sat

    def check_linear_conflict(self, witness):
        conflict = True
        for i in range(n_vars):
            if not witness[i].is_integer():
                aux = witness[i]
                for w in [0,1]:
                    witness[i] = w
                    if (self.lp_solver.A_ub @ witness <= self.lp_solver.y_ub).all():
                        conflict = False
                        break
                witness[i] = aux

        return conflict

