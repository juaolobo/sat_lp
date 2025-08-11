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

    def solve_linear(self, switch=False, fixing={}):

        self.lp_solver.fixing = fixing
        self.lp_solver.create_lp()

        if switch == True:
            self.lp_solver.switch()

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

    def optimize(self, verbose=False):

        it = 0
        witness, res = self.solve_linear()
        n_vars = self.cnf_handler.n_vars
        m_clauses = self.cnf_handler.m_clauses
        stops = 0
        solution_history = set()
        test = set()
        tracking = 0

        if verbose:
            self.print_verbose(witness)

        c = self.lp_solver.c.copy()
        self.lp_solver.restart(last_coefs=c, last_witness=witness)
        min_res = res
        monotonic = [c]
        
        while not self.lp_solver.verify(witness):
            
            last_witness = witness
            if verbose:
                print("------------------------------------------------------------")
                print(f"NUMBER OF DIFFERENT SOLUTIONS: {len(solution_history)} {len(test)}, NUMBER OF CURRENT CYCLE ITERATIONS: {tracking}")

                self.print_verbose(witness)

            witness, res = self.solve_linear()
            # INFEASBLE
            if witness is None:
                return None

            solution_history.add(tuple(self.lp_solver.c))
            test.add(tuple(witness))
            tracking += 1
            self.linear_it += 1

            if len(solution_history) < tracking:
                sat_idx, unsat_idx = self.lp_solver.get_active_clauses(last_witness)
                dmacs_witness = self.lp_solver.coefs_to_witness(last_witness)

                new_clauses = []
                self.bool_solver.restart()
                for idx in unsat_idx:
                    # pick a variable to satisfy the clause
                    # resolve conflict
                    # learn clauses
                    learned = self.bool_solver.expand_and_learn(dmacs_witness, idx)
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
                test = set()
                tracking = 0
                self.boolean_it += 1

                # c = self.lp_solver.c.copy()
                self.lp_solver.restart()

            else:
                c = self.lp_solver.c.copy()
                self.lp_solver.restart(last_coefs=c, last_witness=witness)
                if res < min_res:
                    monotonic.append(c)
                    min_res = res


        if verbose:
            print("------------------------------------------------------------")
            print(f"NUMBER OF DIFFERENT SOLUTIONS: {len(solution_history)} {len(test)}, NUMBER OF CURRENT CYCLE ITERATIONS: {tracking}")

            self.print_verbose(witness)


        print("------------------------------------------------------------")
        print(f"FINAL NUMBER OF ITERATIONS: {self.linear_it}")

        print(monotonic)

        return witness


    def generate_cut_symm_seq(self, verbose=False):

        n_vars = self.cnf_handler.n_vars
        witness = None # placeholder for entering the loop
        fixing_one = {}
        fixing_zero = {}

        while not self.lp_solver.verify(witness):

            # MAX-0
            # solve opt for maximum number of zeros: min sum(z_i)
            # switch obj function back to MAX-0 and fix ones (and zeros) found in previous step
            # print(f"RUNNING MAX-ZERO WITH FIXING_ONE: {fixing_one}")
            witness_zero, res_zero = self.solve_linear(switch=False, fixing=fixing_one)

            # INFEASIBLE => UNSAT
            if witness_zero is None:
                return None, None

            fixing_zero = self.extract_fixing(witness_zero)
            # MAX-1
            # solve opt for maximum number of ones: min -sum(z_i)
            # switch obj function to MAX-1 and fix zeros found in previous step
            # print(f"RUNNING MAX-ONE WITH FIXING_ZERO: {fixing_zero}")
            witness_one, res_one = self.solve_linear(switch=True, fixing=fixing_zero)

            # INFEASIBLE => UNSAT
            if witness_one is None:
                return None, None

            fixing_one = self.extract_fixing(witness_one)
            # detects degeneracy (no boolean coordinates) => lower optimization dimension to 1

            if fixing_one == fixing_zero:
                # check unsat via randomized projection
                fixing, decision = self.randomized_projection(current_fixing=fixing_zero)
                # UNSAT
                if fixing == None:
                    return None, None

                # unable to improve current cut
                if fixing == {}:
                    return fixing_zero, decision

                fixing_one = fixing

            witness = witness_one

        # SAT
        fixing = self.extract_fixing(witness)
        return fixing, None

    def generate_cut_symm_parallel(self, verbose=False):

        n_vars = self.cnf_handler.n_vars
        witness_zero = None # placeholder for entering the loop
        witness_one = None # placeholder for entering the loop
        fixing_one = {}
        fixing_zero = {}
        zero_sat = False
        one_sat = False 

        # add randomized projection at the start once to encourage exploration
        self.lp_solver.fixing = {}
        possible_picks = [i for i in range(n_vars)]
        decision = np.random.choice(possible_picks, 1)[0].item()
        self.lp_solver.create_lp()

        self.lp_solver.set_coefs_for_projection(decision, positive=True)
        witness, res = self.lp_solver.solve()
        fixing_zero = self.extract_fixing(witness)

        self.lp_solver.set_coefs_for_projection(decision, positive=False)
        witness, res = self.lp_solver.solve()
        fixing_one = self.extract_fixing(witness)

        while zero_sat is False or one_sat is False:
            
            # MAX-0
            # solve opt for maximum number of zeros: min sum(z_i)
            # switch obj function back to MAX-0 and fix ones (and zeros) found in previous step
            print(f"RUNNING MAX-ZERO WITH FIXING_ONE: {fixing_one}")
            witness_zero, res_zero = self.solve_linear(switch=False, fixing=fixing_one)
            # INFEASIBLE => UNSAT
            if witness_zero is None:
                return [(None, None)]

            # MAX-1
            # solve opt for maximum number of ones: min -sum(z_i)
            # switch obj function to MAX-1 and fix zeros found in previous step
            print(f"RUNNING MAX-ONE WITH FIXING_ZERO: {fixing_zero}")
            witness_one, res_one = self.solve_linear(switch=True, fixing=fixing_zero)

            # INFEASIBLE => UNSAT
            if witness_one is None:
                return [(None, None)]

            zero_sat = self.lp_solver.verify(witness_zero)
            one_sat = self.lp_solver.verify(witness_one)

            if zero_sat is True or one_sat is True:
                break

            _fixing_zero = self.extract_fixing(witness_zero)
            _fixing_one = self.extract_fixing(witness_one)

            # detects degeneracy (no new boolean coordinates) => lower optimization dimension to 1
            # produce two cuts or learn a new coordinate and continue the optimization

            if fixing_zero == _fixing_one and fixing_one == _fixing_zero:
                possible_cuts = []
                # _fixing mutates fixing_{zero, one} in this loop
                for _fixing in [fixing_one, fixing_zero]:
                    # check unsat via randomized projection
                    fixing, decision = self.randomized_projection(current_fixing=_fixing)

                    # UNSAT
                    if fixing == None:
                        return [(None, None)]

                    # unable to improve current cut
                    if fixing == {}:
                        possible_cuts.append((_fixing, decision))

                    else:
                        if _fixing == fixing_zero:
                            fixing_zero = fixing
                        else:
                            fixing_one = fixing

                if len(possible_cuts) == 2:
                    return possible_cuts

            else:
                fixing_zero = _fixing_zero
                fixing_one = _fixing_one
            
        # SAT
        witness = witness_one if one_sat == True else witness_zero
        fixing = self.extract_fixing(witness)

        return [(fixing, None)]

    def generate_cut_symm_no_fixing(self, verbose=False):

        n_vars = self.cnf_handler.n_vars
        witness_zero = None # placeholder for entering the loop
        witness_one = None # placeholder for entering the loop
        fixing_one = {}
        fixing_zero = {}
        zero_sat = False
        one_sat = False

        while zero_sat is False or one_sat is False:
            
            # MAX-0
            # solve opt for maximum number of zeros: min sum(z_i)
            # switch obj function back to MAX-0 and fix ones (and zeros) found in previous step
            print(f"RUNNING MAX-ZERO WITH FIXING_ONE: {fixing_one}")
            self.lp_solver.create_lp()
            if witness_zero is not None:
                c = [1 if i != 1 else 0 for i in witness_one]
            else:
                c = np.ones(n_vars)

            self.lp_solver.c = c
            witness_zero, res_zero = self.lp_solver.solve()
            # INFEASIBLE => UNSAT
            if witness_zero is None:
                return None, None

            # MAX-1
            # solve opt for maximum number of ones: min -sum(z_i)
            # switch obj function to MAX-1 and fix zeros found in previous step
            print(f"RUNNING MAX-ONE WITH FIXING_ZERO: {fixing_zero}")
            self.lp_solver.create_lp()
            if witness_one is not None:
                c = [-1 if i != 0 else 0 for i in witness_zero]
            else:
                c = np.ones(n_vars)

            self.lp_solver.c = c
            witness_one, res_one = self.lp_solver.solve()

            print(witness_zero)
            print(witness_one)
            # INFEASIBLE => UNSAT
            if witness_one is None:
                return None, None

            zero_sat = self.lp_solver.verify(witness_zero)
            one_sat = self.lp_solver.verify(witness_one)

            if zero_sat is True or one_sat is True:
                break

            # detects degeneracy (no new boolean coordinates) => lower optimization dimension to 1
            # produce two cuts or learn a new coordinate and continue the optimization
            
        # SAT
        witness = witness_one if one_sat == True else witness_zero
        fixing = self.extract_fixing(witness)

        return (fixing, None)

    def test(self):
        cut, decision = self.generate_cut_symm_no_fixing()

        return self.cut_to_linear(cut)

    def symmetric_opt_seq(self):

        n_vars = self.cnf_handler.n_vars
        while True:
            cut, decision = self.generate_cut_symm_seq()
            if cut is None:
                return None

            if len(cut) == n_vars:
                witness = self.cut_to_linear(cut)
                if self.lp_solver.verify(witness):
                    return witness

            else:
                # learn via conflict
                learned_pos = self.bool_solver.expand_and_learn(self.cut_to_witness(cut), decision)
                learned_neg = self.bool_solver.expand_and_learn(self.cut_to_witness(cut), -decision)
                new_clauses = set(learned_pos + learned_neg)
                print(new_clauses)
                for c in new_clauses:
                    self.cnf_handler.add_clause(c)
                    self.cnf_handler.learnt_clauses += 1

    def symmetric_opt_parallel(self, verbose=False):

        generate_cut = self.generate_cut_symm_parallel
        while True:

            possible_cuts = generate_cut()
            
            if len(possible_cuts) == 2: 
                new_clauses = []
                for cut, decision in possible_cuts:
                    # learn via conflict
                    learned_pos = self.bool_solver.expand_and_learn(self.cut_to_witness(cut), decision)
                    learned_neg = self.bool_solver.expand_and_learn(self.cut_to_witness(cut), -decision)
                    new_clauses += set(learned_pos + learned_neg)

                for c in new_clauses:
                    self.cnf_handler.add_clause(c)
                    self.cnf_handler.learnt_clauses += 1

            else:
                cut, _ = possible_cuts[0]

                if cut is None:
                    return None
                
                witness = self.cut_to_linear(cut)
                if self.lp_solver.verify(witness):
                    return witness


    def generate_feas_cut(self):

        fixing = {}
        n_vars = self.cnf_handler.n_vars
        while True:
            witness, res = self.solve_linear(fixing=fixing)
            _fixing = self.extract_fixing(witness)
            
            if fixing == _fixing:

                if len(_fixing) == n_vars:
                    return fixing, None
                print(_fixing)
                _fixing, decision = self.randomized_projection(current_fixing=_fixing)

                # UNSAT
                if _fixing == None:
                    return None, None

                # unable to improve current cut
                if _fixing == {}:
                    return fixing, decision

            fixing = _fixing

        return fixing

    def feas_opt(self):
        n_vars = self.cnf_handler.n_vars
        while True:
            cut, decision = self.generate_feas_cut()
            if cut is None:
                return None

            if len(cut) == n_vars:
                witness = self.cut_to_linear(cut)
                if self.lp_solver.verify(witness):
                    return witness

            else:
                # learn via conflict
                learned_pos = self.bool_solver.expand_and_learn(self.cut_to_witness(cut), decision)
                learned_neg = self.bool_solver.expand_and_learn(self.cut_to_witness(cut), -decision)
                new_clauses = set(learned_pos + learned_neg)
                print(new_clauses)
                for c in new_clauses:
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

    def extract_fixing(self, witness):
        n_vars = self.cnf_handler.n_vars
        fixing = {
            i+1: witness[i].item()
            for i in range(n_vars) 
            if witness[i].is_integer()
        }

        return fixing

    def randomized_projection(self, current_fixing={}):


        # aka randomized projection via optimization
        n_vars = self.cnf_handler.n_vars
        self.lp_solver.fixing = current_fixing
        self.lp_solver.create_lp()

        # pick a variable from the unassigned ones to try expansion
        idxs = [i for i in range(n_vars) if i+1 not in current_fixing.keys()]
        decision = np.random.choice(idxs, 1)[0]

        self.lp_solver.set_coefs_for_projection(decision, positive=True)
        witness, res = self.lp_solver.solve()
        fixing = self.extract_fixing(witness)

        # if positive optimization did not lead to a conflict
        if len(fixing) > len(current_fixing):
            return fixing, -1
                
        self.lp_solver.set_coefs_for_projection(decision, positive=False)
        witness, res = self.lp_solver.solve()
        fixing = self.extract_fixing(witness)

        # if positive optimization led to a conflict and negative didnt
        if len(fixing) > len(current_fixing):
            return fixing, -1

        # if positive and negative led to a conflict
        # and the initial fixing is empty => formula is unsat
        if current_fixing == {}:
            return None, None

        # if both decisions led to a conflict, but the fixing was not empty, 
        # we cant refine the current cut into a boolean solution
        # we return the decision to learn via the conflict

        return {}, abs(decision)+1

    def check_unsat(self, current_fixing={}):

        n_vars = self.cnf_handler.n_vars
        witness = [xi if current_fixing[xi] == 1 else -xi for xi in current_fixing.keys()]
        idxs = [i for i in range(n_vars) if i+1 not in current_fixing.keys()]

        # pick a variable from the unassigned ones to try expansion
        decision = np.random.choice(idxs, 1)[0] + 1

        self.bool_solver.restart()
        expanded_pos = self.bool_solver.expand(witness, decision)

        # if positive decision did not lead to a conflict
        if expanded_pos is not None:
            return current_fixing, decision

        self.bool_solver.restart()
        expanded_neg = self.bool_solver.expand(witness, -decision)

        # if positive decision decision led to a conflict and negative didnt
        if expanded_neg is not None:
            return current_fixing, -decision

        # if both decisions led a conflict, but no variable was assigned prior
        if current_fixing == {}:
            return None, None

        # if both decisions led to a conflict, but the fixing was not empty, 
        # we cant refine the current cut into a boolean solution
        # we return the decision to learn via the conflict
        return {}, decision

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
