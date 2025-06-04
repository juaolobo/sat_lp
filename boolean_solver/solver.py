import numpy as np
import random
import time
from cdcl.clause import Clause
from cdcl.formula import Formula
from cdcl.implication_graph import ImplicationGraph
from cdcl.dimacs_parser import parse

class Solver: 
    def __init__(self, input_cnf_file, verbose):
        self.assert_mode = False
        self.verbose = verbose
        self.list_clause, self.nvars = parse(input_cnf_file, self.verbose)
        self.formula = Formula(self.list_clause)
        self.graph = ImplicationGraph()
        self.decision_level = 0
        self.nb_clauses = len(self.list_clause)
        self.nb_learnt_clause = 0
        self.nb_decisions = 0
        self.restart_count = 0
        self.conflict_count = 0
        self.analysis_count = 0
        self.restart_rate = 100
        self.is_sat = 0 
        self.conflict = None
    
    def restart(self):
        # TODO : Implement other restart mechanisms
        self.formula = Formula(self.list_clause)
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
                pool_literal = w.clause
                i = 0

        return w, w.get_backtrack_level()

    def pick_branching_variable(self):
        
        # if len(self.hypotheses) > 0:
        #     idx = np.random.choice(range(len(self.hypotheses)), replace=False)
        #     decision = self.hypotheses.pop(idx)
        #     if (decision in self.graph.assigned_vars) or (-decision in self.graph.assigned_vars):
        #         return self.pick_branching_variable()
        #     else:
        #         return decision

        # else: 
        ## Most frequent var first
        # counter = self.formula.get_counter()
        counter = {}
        for clause in self.formula.formula:
            for literal in clause.clause[:clause.size]:
                if literal in counter:
                    counter[literal] += 1
                else:
                    counter[literal] = 1 
        assert len(counter) > 0

        pool_literal = list(counter.keys())
        # pool_literal = sorted(list(counter.keys()))
        decision = -1
        i = 0
        while i < len(pool_literal):
            decision = pool_literal[i]
            if decision not in self.graph.assigned_vars and -decision not in self.graph.assigned_vars:
                break
            i += 1

        assert decision not in self.graph.assigned_vars
        assert -decision not in self.graph.assigned_vars
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
        
        witness = None
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
        

        return self.formula.get_value, witness    
            

    def propagate_linear(self, linear_sol):

        _linear_sol = linear_sol
        self.is_sat, self.conflict =  self.formula.unit_propagate(self.decision_level, self.graph)
        for lit in _linear_sol:
            self.graph.add_node(lit, None, 0)
            self.is_sat, self.conflict = self.formula.bcp(lit, 0, self.graph)
            if self.is_sat != 0:
                # leave this here if it ever happens
                print("SUPRISE MOTHERFUCKER")
                break

        self.is_sat, self.conflict = self.formula.unit_propagate(1, self.graph)
        self.decision_level += 1
        conflict = None
        learnt_clause = None
        m_simpl = -1
        m_clauses = len(self.formula.formula)
        while m_simpl < m_clauses:

            # simplify all possible clauses by deduction ([~p, q], [p, q] -> q)
            m_simpl = 0
            for c in self.formula.formula:
                if c.size == 2:
                    inf, clause = self.formula.simplify_pair(c, self.graph, self.decision_level)
                    if inf is not None and clause is not None:
                        self.graph.add_node(inf, clause, self.decision_level)
                        self.is_sat, self.conflict = self.formula.bcp(inf, self.decision_level, self.graph)
                        if self.is_sat == -1:
                            break
                m_simpl += 1

            if self.is_sat == 0:
                self.is_sat, self.conflict = self.formula.unit_propagate(self.decision_level, self.graph)

            # solve conflict
            while self.is_sat == -1:
                learnt_clause, backtrack_level = self.conflict_analysis(self.conflict)
                self.formula.add_clause(learnt_clause)
                self.graph.backtrack(backtrack_level)
                self.formula.backtrack(backtrack_level, self.graph)
                for xi in learnt_clause.clause:
                    self.formula.repair(-xi, self.graph)
                    self.graph.remove_node(-xi)

                _linear_sol = [xi for xi in _linear_sol if -xi not in learnt_clause.clause]

                self.is_sat, self.conflict = self.formula.unit_propagate(self.decision_level, self.graph)

        # return once there are no deductions to be made and all conflicts have been resolved
        return self.graph.assigned_vars, self.formula

    def extend_solution(self):
        print("WELL")
        
        # assign variables to advance further in the problem (and generate conflicts)
        while self.is_sat == 0:
            
            decision = self.pick_branching_variable()
            self.decision_level += 1
            self.graph.add_node(decision, None, self.decision_level)
            self.is_sat, self.conflict = self.formula.bcp(decision, self.decision_level, self.graph)

        # solve all conflicts generated
        while self.is_sat == -1:

            learnt_clause, backtrack_level = self.conflict_analysis(self.conflict)
            self.formula.add_clause(learnt_clause)
            self.graph.backtrack(backtrack_level)
            self.formula.backtrack(backtrack_level, self.graph)

            self.is_sat, self.conflict = self.formula.unit_propagate(self.decision_level, self.graph)

        return self.graph.assigned_vars, self.formula

"""
-- outline of the algorithm
- get linear solver partial answer
- simplify all clauses
    - fix conflicts if the appear and return to linear solver
    - if there are no conflicts, assign variables until there are conflicts
        - fix conflicts and return to linear solver
"""