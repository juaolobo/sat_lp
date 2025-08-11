import numpy as np
import time
from satlp.boolean_solver import Clause, Formula, ImplicationGraph
from satlp.cnf_loader import CNFLoader

class BooleanSolver: 
    def __init__(self, filename, verbose, cnf_handler=None):
        self.verbose = verbose
        self.cnf_handler = cnf_handler if cnf_handler is not None else CNFLoader(filename)
        self.nvars = self.cnf_handler.n_vars
        self.formula = Formula(self.cnf_handler.clauses)
        self.graph = ImplicationGraph()
        self.decision_level = 0
        self.nb_clauses = len(self.cnf_handler.clauses)
        self.nb_learnt_clause = 0
        self.nb_decisions = 0
        self.restart_count = 0
        self.conflict_count = 0
        self.analysis_count = 0
        self.restart_rate = 100
        self.is_sat = 0 
        self.conflict = None

    def witness_to_linear(self, witness):
        witness = sorted(witness, key=abs)
        solution = [1 if xi > 0 else 0 for xi in witness]

        return solution

    def restart(self):
        self.formula = Formula(self.cnf_handler.clauses)
        self.graph = ImplicationGraph()
        self.decision_level = 0
        self.restart_count += 1
        self.conflict_count = 0
        self.is_sat = 0
        self.conflict = None
            
    def conflict_analysis(self, conflict_clause): 
        # RETURN : learnt clause, and backtrack_level
        w = conflict_clause
        pool_literal = w.clause
        i = 0
        while i < len(pool_literal):
            conflict_literal = pool_literal[i]
            antecedent = self.graph.get_antecedent(-conflict_literal)
            i += 1

            if antecedent is not None:
                w = w.resolution_operate(antecedent, conflict_literal)
                if w is None:
                    return None, -1

                pool_literal = w.clause
                i = 0

        return w, w.get_backtrack_level()
    
    def pick_branching_variable(self):

        ## Most frequent var first
        counter = {}
        for clause in self.formula.formula:
            for literal in clause.clause[:clause.size]:
                if literal in counter:
                    counter[literal] += 1
                else:
                    counter[literal] = 1 
        assert len(counter) > 0

        pool_literal = list(counter.keys())
        decision = -1
        i = 0
        while i < len(pool_literal):
            decision = pool_literal[i]
            if decision not in self.graph.assigned_vars and -decision not in self.graph.assigned_vars:
                break
            i += 1

        assert decision not in self.graph.assigned_vars
        assert -decision not in self.graph.assigned_vars
        # decision = -decision if np.random.uniform(0, 1) >= 0.5 else decision

        return decision

    def is_all_assigned(self):
        return self.nvars == len(self.graph.assigned_vars)

    def solve(self): 
        stop = False
        initial_time = time.time()
        self.is_sat, self.conflict =  self.formula.unit_propagate(self.decision_level, self.graph)
        if self.verbose:
            print('=====================[  Search Statistics ]=====================')
            
        while self.is_sat == 0 and not stop: 
            assert self.formula.get_value() == self.is_sat
            assert self.conflict is None
            assert not self.is_all_assigned() 

            decision = self.pick_branching_variable()
            self.nb_decisions += 1
            self.decision_level += 1
            self.graph.add_node(decision, None, self.decision_level)
            self.is_sat, self.conflict = self.formula.bcp(decision, self.decision_level, self.graph)

            if self.is_sat == 0:
                assert not self.is_all_assigned()
                self.is_sat, self.conflict = self.formula.unit_propagate(self.decision_level, self.graph)

            if self.is_sat == 0:
                assert not self.is_all_assigned()

            if self.is_sat == 1:
                assert self.formula.get_value() == self.is_sat
                break

            while self.is_sat == -1 and not stop:
                assert self.conflict is not None
                learnt_clause, backtrack_level = self.conflict_analysis(self.conflict)

                if backtrack_level == -100:
                    self.restart()      
                elif backtrack_level == -1 :
                    self.is_sat = -1 
                    stop = True
                else:
                    self.formula.add_clause(learnt_clause)
                    self.cnf_handler.add_clause(learnt_clause.clause)
                    self.nb_learnt_clause += 1
                    self.graph.backtrack(backtrack_level)
                    self.formula.backtrack(backtrack_level, self.graph)
                    self.decision_level = backtrack_level
                    self.is_sat, self.conflict = self.formula.unit_propagate(self.decision_level, self.graph)
                    if self.is_sat == 0: 
                        assert not self.is_all_assigned()

                ## If too much conflicts, RESTART IT NOW !
                self.conflict_count += 1
                if self.conflict_count > self.restart_rate:
                    self.restart()

        assert self.is_sat != 0
        assert self.is_sat == self.formula.get_value()
        if self.is_all_assigned():
            print('All vars assigned !')
        else: 
            print('Early quit !')

        print('Restart: ', self.restart_count)
        print('Learnt clauses: ', self.nb_learnt_clause)
        print('Decisions: ', self.nb_decisions)
        print('CPU time: {0:10.6f}s '.format(time.time()-initial_time))
        
        witness = self.graph.assigned_vars

        if stop: 
            assert self.is_sat == -1
            print('UNSAT')
        elif not stop and self.is_sat == 1:
            stop = self.graph.assigned_vars
            print('SAT')
        elif not stop and self.is_sat == -1:
            print('UNSAT')
        else: #But practically, this should not happen !
            print('UNRESOLVED !')
        
        # solution was found without assigning all variable
        if len(witness) < self.nvars:
            for xi in range(1, self.nvars+1):
                if xi not in witness and -xi not in witness:
                    witness.append(xi)

        return witness    

    def propagate_linear(self, linear_sol):

        self.is_sat, self.conflict =  self.formula.unit_propagate(self.decision_level, self.graph)
        self.fix_variables(linear_sol)
        self.is_sat, self.conflict = self.formula.unit_propagate(1, self.graph)
        self.decision_level += 1
        conflict = None
        learnt_clause = None
        m_simpl = -1
        m_clauses = len(self.formula.formula)
        while m_simpl < m_clauses:

            # simplify all possible clauses by deduction ([~p, q], [p, q] -> q)
            m_simpl = 0
            for i in range(m_clauses):
                c = self.formula.formula[i]
                if c.size == 2:
                    inf, clause = self.formula.simplify_pair(i, self.graph, self.decision_level)
                    if inf is not None and clause is not None:
                        self.graph.add_node(inf, clause, self.decision_level)
                        self.is_sat, self.conflict = self.formula.bcp(inf, self.decision_level, self.graph)
                        if self.is_sat == -1:
                            break
                m_simpl += 1

            if self.is_sat == 0:
                self.is_sat, self.conflict = self.formula.unit_propagate(self.decision_level, self.graph)

            new_clauses = []
            # solve conflict
            while self.is_sat == -1:

                learnt_clause, backtrack_level = self.conflict_analysis(self.conflict)

                # detected unsolvable conflict => UNSAT
                if learnt_clause is None:
                    return None, None
                

                self.formula.add_clause(learnt_clause)
                self.graph.backtrack(backtrack_level)
                self.formula.backtrack(backtrack_level, self.graph)

                for xi in learnt_clause.clause:
                    self.formula.repair(-xi, self.graph)
                    self.graph.remove_node(-xi)

                new_clauses.append(learnt_clause.clause)
                self.is_sat, self.conflict = self.formula.unit_propagate(self.decision_level, self.graph)

        # return once there are no deductions to be made and all conflicts have been resolved
        return self.graph.assigned_vars, new_clauses

    def extend_solution(self):
        
        new_clauses = []
        # assign variables to advance further in the problem (and generate conflicts)
        while self.is_sat == 0:
            
            decision = self.pick_branching_variable()
            self.decision_level += 1
            self.graph.add_node(decision, None, self.decision_level)
            self.is_sat, self.conflict = self.formula.bcp(decision, self.decision_level, self.graph)

        # solve all conflicts generated
        while self.is_sat == -1:

            learnt_clause, backtrack_level = self.conflict_analysis(self.conflict)
            # detected unsolvable conflict => UNSAT
            if learnt_clause is None:
                return None, None

            self.formula.add_clause(learnt_clause)
            self.graph.backtrack(backtrack_level)
            self.formula.backtrack(backtrack_level, self.graph)
            new_clauses.append(learnt_clause.clause)

            self.is_sat, self.conflict = self.formula.unit_propagate(self.decision_level, self.graph)

        return self.graph.assigned_vars, new_clauses

    def pick_sat_var(self, clause):
        lit_pool = clause.clause[:clause.size]
        
        return lit_pool
        decision = np.random.choice(lit_pool, 1).item()

        assert decision not in self.graph.assigned_vars
        return decision

    def fix_variables(self, linear_sol):
        
        for lit in linear_sol:
            if lit not in self.graph.assigned_vars:
                self.graph.add_node(lit, None, 0)
                self.is_sat, self.conflict = self.formula.bcp(lit, 0, self.graph)

                if self.is_sat != 0:
                    # leave this here if it ever happens
                    print("SUPRISE MOTHERFUCKER", self.is_sat)
                    # self.conflict.print_info()
                    return

    def _expand_and_learn(self, witness, idx):
        
        clause = self.formula.formula[idx]
        # current bug: fixed variables are not remaining fixed in later iterations
        self.fix_variables(witness)

        # should break produce a conflict
        new_clauses = []
        decisions = self.pick_sat_var(clause)

        for decision in decisions:
            self.decision_level += 1

            self.graph.add_node(decision, None, self.decision_level)
            self.is_sat, self.conflict = self.formula.bcp(decision, self.decision_level, self.graph)

            if self.is_sat == 0:
                self.is_sat, self.conflict = self.formula.unit_propagate(self.decision_level, self.graph)

            # solve all conflicts generated
            while self.is_sat == -1:
                learnt_clause, backtrack_level = self.conflict_analysis(self.conflict)
                # detected unsolvable conflict => UNSAT
                if learnt_clause is None:
                    breakpoint()
                    return None

                new_clauses.append(learnt_clause.clause)
                self.formula.add_clause(learnt_clause)
                self.graph.backtrack(backtrack_level)
                self.formula.backtrack(backtrack_level, self.graph)
                    
                self.is_sat, self.conflict = self.formula.unit_propagate(self.decision_level, self.graph)

            self.restart()
            self.fix_variables(witness)

        self.restart()
        return new_clauses

    def expand_and_learn(self, witness, decision):
        
        self.restart()
        self.fix_variables(witness)

        self.decision_level += 1
        self.graph.add_node(decision, None, self.decision_level)
        self.is_sat, self.conflict = self.formula.bcp(decision, self.decision_level, self.graph)

        if self.is_sat == 0:
            self.is_sat, self.conflict = self.formula.unit_propagate(self.decision_level, self.graph)

        # solve all conflicts generated
        new_clauses = []
        while self.is_sat == -1:
            learnt_clause, backtrack_level = self.conflict_analysis(self.conflict)
            # detected unsolvable conflict => UNSAT
            if learnt_clause is None:
                return None

            new_clauses.append(tuple(learnt_clause.clause))
            self.formula.add_clause(learnt_clause)
            self.graph.backtrack(backtrack_level)
            self.formula.backtrack(backtrack_level, self.graph)
                
            self.is_sat, self.conflict = self.formula.unit_propagate(self.decision_level, self.graph)

        return new_clauses


"""
-- outline of the algorithm
- get linear solver partial answer
- simplify all clauses
    - fix conflicts if the appear and return to linear solver
    - if there are no conflicts, assign variables until there are conflicts
        - fix conflicts and return to linear solver
"""