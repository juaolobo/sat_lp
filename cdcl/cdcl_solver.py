from cdcl import parse
import numpy as np
import random
import time
from cdcl.cnf_data_structure import Clause, CNF_Formula, Implication_Graph

class CDCL_Solver: 
    def __init__(self, input_cnf_file, verbose, hypotheses=[]):
        self.assert_mode = False
        self.verbose = verbose
        self.list_clause, self.nvars = parse(input_cnf_file, self.verbose)
        self.formula = CNF_Formula(self.list_clause)
        self.graph = Implication_Graph()
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
        self.hypotheses = hypotheses
    
    def restart(self):
        # TODO : Implement other restart mechanisms
        self.formula = CNF_Formula(self.list_clause)
        self.graph = Implication_Graph()
        self.decision_level = 0
        self.restart_count += 1
        self.conflict_count = 0
        self.is_sat = 0
        self.conflict = None
    
    def conflict_analysis(self, conflict_clause): 
        # RETURN : learnt clause, and backtrack_level
        self.analysis_count += 1
        if self.analysis_count >= 100:
            # If conflict analysis is called too long => restart 
            self.analysis_count = 0
            return None, -100
        w = conflict_clause
        pool_literal = w.literal_at_level(self.decision_level)
        assert len(pool_literal) > 0
        if len(pool_literal) == 1:
            self.analysis_count = 0
            return w, w.get_backtrack_level()
        else: 
            i = 0
            antecedent = None
            while i < len(pool_literal) and antecedent is None:
                conflict_literal = pool_literal[i]
                antecedent = self.graph.get_antecedent(-conflict_literal)
                i += 1
            if antecedent is not None:
                assert antecedent is not None
                w = w.resolution_operate(antecedent, conflict_literal)
                return self.conflict_analysis(w)
            
            else:
                # should happen only for hypothesis
                return w, None

    def pick_branching_variable(self):

        # if len(self.hypothesis) > 0:
        #     idx = np.random.choice(range(len(self.hypothesis)), replace=False)
        #     decision = self.hypothesis.pop(idx)
        #     if (decision in self.graph.assigned_vars) or (-decision in self.graph.assigned_vars):
        #         return self.pick_branching_variable()

        # else:
        ## Most frequent var first
        counter = self.formula.get_counter()
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
        
        if stop: 
            assert self.is_sat == -1
            print('UNSAT')
        elif not stop and self.is_sat == 1:
            print('SAT')
        elif not stop and self.is_sat == -1:
            print('UNSAT')
        else: #But practically, this should not happen !
            print('UNRESOLVED !')

        ## Check it
        assigned_vars = self.graph.assigned_vars
        test_formula = CNF_Formula(self.list_clause)
        test_formula.formula = [Clause(c) for c in self.list_clause if len(c) > 0]
        test_graph = Implication_Graph()
        for i, literal in enumerate(assigned_vars):
            test_formula.bcp(literal, i, test_graph)

        if self.is_sat == 1:
            assert test_formula.get_value() == 1

        print('Verified by (re)propagating assignments on original clauses !')
            
        return self.formula.get_value    

    def solve_for_conflict(self):
        stop = False
        decision_made = False
        self.is_sat, self.conflict =  self.formula.unit_propagate(self.decision_level, self.graph)
        # change to force solver to make a choice
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
                # if backtrack_level == -100:
                #     self.restart()      
                if backtrack_level == -1 or backtrack_level == -100:
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

                if self.is_sat != -1:
                    return self.graph.assigned_vars, self.formula

        return None, None   

    def solve_for_conflict2(self):
        stop = False
        decision_made = False
        for lit in self.hypotheses:
            self.graph.add_node(lit, None, 0)
            self.is_sat, self.conflict = self.formula.bcp(lit, self.decision_level, self.graph)
            if self.is_sat != 0:
                # leave this here if it ever happens
                print("SUPRISE MOTHERFUCKER")
                break         
        
        breakpoint()
        if self.is_sat == 0:
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

        while self.is_sat == -1 and not stop:
            assert self.conflict is not None
            learnt_clause, backtrack_level = self.conflict_analysis(self.conflict)
            # if backtrack_level == -100:
            #     self.restart()      

            if backtrack_level == 0:
                # then it is a hypothesis
                conflict_lits = self.conflict.literal_at_level(0)
                exit()
                # self.conflict.size = len(conflict_lits)
                # self.conflict.decision_level = [xi if xi > 0 else -1 for xi in self.conflict.decision_level]
                # self.conflict.value = 0
                # for lit in conflict_lits:
                #     self.graph.remove_node(lit)

                # self.conflict.update(self.graph)
                # self.conflict.pick_new_ref()
                # self.is_sat, self.conflict = self.formula.unit_propagate(self.decision_level, self.graph)

                # if self.conflict is None:
                #     return self.graph.assigned_vars, self.formula

            elif backtrack_level == -1 or backtrack_level == -100:
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

            if self.is_sat != -1:
                return self.graph.assigned_vars, self.formula

        return None, None
            

