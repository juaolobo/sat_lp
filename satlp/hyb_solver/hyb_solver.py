from satlp.linear_solver import (
    SATasLPOptimization,
    SATasLPFeasibility,
    SATasLPOptimizationDual,
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
    ):
        self.filename = filename
        self.fixing = {}
        self.cnf_handler = CNFLoader(filename)
        self.lp_solver = lp_solver(fixing=self.fixing, filename=filename, method=method, cnf_handler=self.cnf_handler)
        self.bool_solver = BooleanSolver(filename, verbose=0, cnf_handler=self.cnf_handler)
        self.history = []
        self.linear_it = 0
        self.boolean_it = 0
        self.wp_it = 0
        

    def print_verbose(self, witness):
        n_vars = self.cnf_handler.n_vars
        dmacs = self.lp_solver.coefs_to_witness(witness)
        print(f"X: {witness[:n_vars]}")
        print(f"X_+: {witness[n_vars:2*n_vars]}")
        print(f"C_+: {self.lp_solver.c[n_vars:2*n_vars]}")
        print(f"X_-: {witness[2*n_vars:]}")
        print(f"C_-: {self.lp_solver.c[2*n_vars:]}")

        print(f"WITNESS: {sorted(list(dmacs), key=abs)}")

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

    def solve_linear(self, switch=False, fixing={}, feas=False):

        self.lp_solver.fixing = fixing
        self.lp_solver.create_lp()

        if switch == True:
            self.lp_solver.switch()

        if feas:
            n_vars = self.cnf_handler.n_vars
            self.lp_solver.c = np.zeros(n_vars)


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

    def generate_cut_symm(self, verbose=False):

        n_vars = self.cnf_handler.n_vars
        witness = None # placeholder for entering the loop
        fixing_one = {}
        fixing_zero = {}
        print("------------------------------------------------------------------------------")
        conflicts = []
        while True:

            # MAX-0
            # solve opt for maximum number of zeros: min sum(z_i)
            # switch obj function back to MAX-0 and fix ones (and zeros) found in previous step
            print(f"RUNNING MAX-ZERO WITH FIXING_ONE: {fixing_one}")
            switch = True if np.random.uniform() > 0.5 else False
            witness_zero, res_zero = self.solve_linear(switch=switch, fixing=fixing_one)
            self.linear_it += 1
            # INFEASIBLE => UNSAT
            if witness_zero is None:
                return None, None

            if self.lp_solver.verify(witness_zero):
                witness = witness_zero
                break

            fixing_zero = self.extract_fixing(witness_zero)

            # MAX-1
            # solve opt for maximum number of ones: min -sum(z_i)
            # switch obj function to MAX-1 and fix zeros found in previous step
            print(f"RUNNING MAX-ONE WITH FIXING_ZERO: {fixing_zero}")
            switch = True if np.random.uniform() > 0.5 else False
            witness_one, res_one = self.solve_linear(switch=(switch), fixing=fixing_zero)
            self.linear_it += 1

            # INFEASIBLE => UNSAT
            if witness_one is None:
                return None, None

            if self.lp_solver.verify(witness_one):
                witness = witness_one
                break

            fixing_one = self.extract_fixing(witness_one)
            # detects degeneracy (no boolean coordinates) => lower optimization dimension to 1

            if fixing_one == fixing_zero:
                # check unsat via weak branching
                witness, decision = self.weak_projection(current_fixing=fixing_zero)
                self.wp_it += 1
                fixing = self.extract_fixing(witness)

                if abs(decision) not in fixing.keys():
                    
                    conflicts.append((fixing, decision))

                    _witness, _decision = self.weak_projection(
                        current_fixing=fixing_zero, 
                        decision=-decision,
                    )
                    self.wp_it += 1
                    _fixing = self.extract_fixing(_witness)

                    if len(_fixing) == n_vars:
                        return _fixing, []

                    # UNSAT or conflict
                    if abs(_decision) not in _fixing.keys():
                        return _fixing, conflicts

                    fixing = _fixing

                fixing_one = fixing

            witness = witness_one


        # SAT
        fixing = self.extract_fixing(witness)
        return fixing, []

    def generate_cut_via_weak_projection(self, verbose=False, init_fixing={}):

        n_vars = self.cnf_handler.n_vars
        fixing = init_fixing
        print("------------------------------------------------------------------------------")
        print(f"BOOTSTRAPING WITH {init_fixing}")
        conflicts = []
        # np.random.seed(41)
        while True:

            _witness, _decision = self.weak_projection(current_fixing=fixing)
            self.wp_it += 1
            _fixing = self.extract_fixing(_witness)
            print(f"TRIED WEAK BRANCHING ON COORDINATE {_decision}")

            print(f"RESULT: {_fixing}")
            if len(_fixing) == n_vars:
                return _fixing, []


            if abs(_decision) not in _fixing.keys():

                conflicts.append((_fixing,_decision))

                _witness, _decision = self.weak_projection(
                    current_fixing=_fixing, 
                    decision=-_decision
                )
                self.wp_it += 1
                _fixing = self.extract_fixing(_witness)

                print(f"TRIED WEAK BRANCHING ON COORDINATE {_decision}")

                print(f"RESULT: {_fixing}")

                if len(_fixing) == n_vars:
                    return _fixing, []
                # UNSAT or conflict
                if abs(_decision) not in _fixing.keys():
                    return _fixing, conflicts

            fixing = _fixing

        # SAT
        fixing = self.extract_fixing(witness)
        return fixing, []

    def generate_feas_cut(self):

        fixing = {}
        n_vars = self.cnf_handler.n_vars
        conflicts = []
        print("------------------------------------------------------------------------------")

        while True:

            print(f"RUNNING FEASIBILITY WITH FIXING: {fixing}")
            witness, res = self.solve_linear(fixing=fixing, feas=True)
            self.linear_it += 1
            _fixing = self.extract_fixing(witness)

            if len(_fixing) == n_vars:
                return _fixing, []
                

            if fixing == _fixing:
                _witness1, _decision = self.weak_projection(current_fixing=_fixing)
                self.wp_it += 1
                _fixing = self.extract_fixing(_witness1)

                if len(_fixing) == n_vars:
                    return _fixing, []

                if abs(_decision) not in _fixing.keys():
                    
                    conflicts.append((_fixing, _decision))

                    print(f"RUNNING WEAK BRANCHING ON COORDINATE {abs(_decision)} AND DIRECTION {-_decision}")
                    _witness, _decision = self.weak_projection(
                        current_fixing=_fixing, 
                        decision=-_decision
                    )
                    self.wp_it += 1

                    _fixing = self.extract_fixing(_witness)
                    print(f"RESULT: {_fixing}")

                    if len(_fixing) == n_vars:
                        return _fixing, []

                    # UNSAT or conflict
                    if abs(_decision) not in _fixing.keys():

                        return _fixing, conflicts

            fixing = _fixing

        return fixing

    def optimize(self, generate_cut, track_history=False):

        n_vars = self.cnf_handler.n_vars
        fixing = {}
        while True:

            cut, conflicts = generate_cut()
            if track_history is True:
                self.history.append(cut)

            # UNSAT
            if len(cut) == 0:
                return None

            if len(cut) == n_vars:
                witness = self.cut_to_linear(cut)
                if self.lp_solver.verify(witness):
                    return witness

            else:
                # learn via conflict
                for _cut, decision in conflicts:

                    learned = self.bool_solver.expand_and_learn(self.cut_to_witness(_cut), decision)
                    if len(learned) == 0:
                        with open("debug.log", "a") as db:
                            db.write(self.filename)
                            db.write('\n')
                            db.write(str(self.cnf_handler.clauses))
                            db.write('\n')
                            db.write(str(conflicts))
                        breakpoint()
                    # UNSAT
                    if learned is None:
                        return None

                    for c in learned:
                        self.cnf_handler.add_clause(c)
                        self.cnf_handler.learnt_clauses += 1
                    

    def cut_to_linear(self, cut):
        n_vars = self.cnf_handler.n_vars
        linear = [cut[i+1] for i in range(n_vars)]

        return linear

    def linear_to_witness(self, linear_sol):
        n_vars = self.cnf_handler.n_vars
        witness = [i+1 if linear_sol[i]== 1 else -i-1 for i in range(n_vars) if linear_sol[i].is_integer()]

        return witness

    def cut_to_witness(self, cut):
        witness = [xi if cut[xi] == 1 else -xi for xi in cut.keys()]
        return witness

    def is_boolean(self, x, eps=1e-8):
        return np.abs(x) < eps or np.abs(1-x) < eps

    def extract_fixing(self, witness):
        n_vars = self.cnf_handler.n_vars
        fixing = {
            i+1: round(witness[i].item())
            for i in range(n_vars) 
            if self.is_boolean(witness[i])
        }

        return fixing

    def weak_projection(self, current_fixing={}, decision=None):

        # aka randomized projection via optimization
        n_vars = self.cnf_handler.n_vars
        self.lp_solver.fixing = current_fixing
        self.lp_solver.create_lp()

        exclude = [xi if current_fixing[xi] == 1 else -xi for xi in current_fixing.keys()]
        if decision is None:
            decision = self.bool_solver.pick_branching_variable(exclude=exclude)

        b = 0 if decision < 0 else 1
        idx = abs(decision)-1

        self.lp_solver.set_coefs_for_projection(idx, direction=b)
        witness, res = self.lp_solver.solve()

        return witness, decision

    def verify(self, witness):
        
        sat = False
        n_vars = self.cnf_handler.n_vars
        if witness is not None:
            sat = self.lp_solver.verify(witness[:n_vars])
            if sat is True:
                print("SATISFIABLE")
                witness = self.linear_to_witness(witness)
                print(f"WITNESS: {witness[:n_vars]}")

            else:
                print("SOMETHING IS WRONG")

            return sat

        print("UNSATISFIABLE")
        return sat
