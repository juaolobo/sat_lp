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
        self.formula = CNF_Formula(self.list_clause, hypotheses)
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
    
    def restart(self, hypotheses=[]):
        # TODO : Implement other restart mechanisms
        self.hypotheses = hypotheses
        self.formula = CNF_Formula(self.list_clause, self.hypotheses)
        self.graph = Implication_Graph()
        
        self.decision_level = 0
        self.restart_count += 1
        self.conflict_count = 0
        self.is_sat = 0
        self.conflict = None
    
    def conflict_analysis(self, conflict_clause): 
        # RETURN : learnt clause, and backtrack_level
        self.analysis_count += 1
        # if self.analysis_count >= 100:
        #     # If conflict analysis is called too long => restart 
        #     self.analysis_count = 0
        #     return None, -100
        w = conflict_clause
        pool_literal = w.literal_at_level(self.decision_level)
        # assert len(pool_literal) > 0
        # if len(pool_literal) == 1:
        if len(pool_literal) == 1:
            conflict_literal = pool_literal[0]
            antecedent = self.graph.get_antecedent(-conflict_literal)
            if antecedent is not None:
                w = w.resolution_operate(antecedent, conflict_literal)

            self.analysis_count = 0
            return w, w.get_backtrack_level()
        else: 
            i = 0
            antecedent = None
            while i < len(pool_literal) and antecedent is None:
                conflict_literal = pool_literal[i]
                antecedent = self.graph.get_antecedent(-conflict_literal)
                i += 1

            assert antecedent is not None
            w = w.resolution_operate(antecedent, conflict_literal)

            return self.conflict_analysis(w)
            
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

        return w            


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
        counter = self.formula.get_counter()
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
        test_formula = CNF_Formula(self.list_clause, hypotheses=[])
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
        for lit in self.hypotheses:
            self.graph.add_node(lit, None, 0)
            if self.is_sat != 0:
                # leave this here if it ever happens
                print("SUPRISE MOTHERFUCKER")
                break

        for lit in self.hypotheses:
            self.is_sat, self.conflict = self.formula.bcp(lit, 0, self.graph)

        self.is_sat, self.conflict = self.formula.unit_propagate(self.decision_level, self.graph)

        decision = self.pick_branching_variable()
        decision = 7 # breaks
        self.decision_level += 1
        self.nb_decisions += 1
        self.graph.add_node(decision, None, self.decision_level)
        self.is_sat, self.conflict = self.formula.bcp(decision, self.decision_level, self.graph)

        breakpoint()
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
                self.graph.backtrack(backtrack_level)
                self.formula.backtrack(backtrack_level, self.graph)
                self.decision_level = backtrack_level
                self.is_sat, self.conflict = self.formula.unit_propagate(self.decision_level, self.graph)

                self.formula.add_clause(learnt_clause)
                new_hypotheses = []
                for xi in self.hypotheses:
                    if (-xi in learnt_clause.clause) or (xi in learnt_clause.clause):
                        self.graph.remove_node(xi)
                    else:
                        new_hypotheses.append(xi)

                # new_hypotheses.append(decision)
                self.hypotheses = new_hypotheses

                # resolved = [xi for xi in self.hypotheses if (-xi not in learnt_clause.clause) and xi not in learnt_clause.clause]

                return self.hypotheses, self.formula
                
            elif backtrack_level == -1 or backtrack_level == -100:
                self.is_sat = -1 
                stop = True

            else:
                print("ELSE")
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

            if self.is_sat == 0:
                # breakpoint()
                return self.hypotheses, self.formula
                # return self.graph.assigned_vars, self.formula

        # self.hypotheses.append(decision)
        return None, None
            

    def solve_for_real(self):

        self.is_sat, self.conflict =  self.formula.unit_propagate(self.decision_level, self.graph)
        for lit in self.hypotheses:
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
            m_simpl = 0
            for c in self.formula.formula:
                
                if c.size == 2:
                    inf, clause = self.formula.simplify_pair(c, self.graph, self.decision_level)
                    if inf is not None and clause is not None:
                        # self.formula.add_clause(clause)
                        self.graph.add_node(inf, clause, self.decision_level)
                        self.is_sat, self.conflict = self.formula.bcp(inf, self.decision_level, self.graph)
                        # self.decision_level += 1
                        if self.is_sat == -1:
                            break
                m_simpl += 1

            if self.is_sat == 0:
                self.is_sat, self.conflict = self.formula.unit_propagate(self.decision_level, self.graph)

            while self.is_sat == -1:
                learnt_clause, backtrack_level = self.conflict_analysis(self.conflict)
                learnt_clause.print_info()
                self.formula.add_clause(learnt_clause)
                self.graph.backtrack(0)
                self.formula.backtrack(0, self.graph)
                # learnt_clause.decision_level = [-1 for _ in range(len(learnt_clause.clause))]
                for xi in learnt_clause.clause:
                    self.formula.repair(-xi, self.graph)
                    self.graph.remove_node(-xi)

                self.hypotheses = [xi for xi in self.hypotheses if -xi not in learnt_clause.clause]

                self.is_sat, self.conflict = self.formula.unit_propagate(self.decision_level, self.graph)

            # bug: some learnt clauses are weird

        return self.graph.assigned_vars, self.formula