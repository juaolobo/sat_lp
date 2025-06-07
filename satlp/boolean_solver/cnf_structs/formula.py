from satlp.boolean_solver.cnf_structs.clause import Clause
import numpy as np

class Formula:
    
    def __init__(self, list_clause):
        self.formula = []
        for c in list_clause:
            if len(c) > 0:
                clause = Clause(c)
                self.formula.append(clause)

        self.value = self.get_value()

    def print_info(self):
        for c in self.formula:
            c.print_info()
        print('[F] Truth value: ', self.value)

    def get_value(self):
        list_values = [c.value for c in self.formula]
        if -1 in list_values:
            return -1
        else:
            if 0 in list_values:
                return 0
            else: 
                return 1 

    def get_counter(self):

        # change this
        counter = {}
        for clause in self.formula:
            for literal in clause.clause[:clause.size]:
                if literal in counter:
                    counter[literal] += 1
                else:
                    counter[literal] = 1 

        ## Lazy counter

        unassigned_refs = []
        for clause in self.formula:
            if clause.size == 1:
                unassigned_refs.append(clause.clause[0])
            elif clause.size > 1:
                for xi in clause.clause[clause.size:]:
                    unassigned_refs.append(xi)

        for literal in unassigned_refs:
            if literal in counter:
                counter[literal] += 1
            else:
                counter[literal] = 1 
        counter = {k: freq for k, freq in sorted(counter.items(), key=lambda item: item[1], reverse= True)}
        assert len(counter) > 0
        return counter

    def is_sat(self):
        return self.value
    
    def simplify_pair(self, pair, graph, decision_level):

        assigned_vars = graph.assigned_vars
        for cl in self.formula:
            if cl.size == 2:
                if cl != pair:
                    _cl = np.array(sorted(cl.clause[:cl.size], key=abs))
                    _pair = np.array(sorted(pair.clause[:pair.size], key=abs))
                    idx = (_cl == _pair)

                    if idx.sum() == 1:
                        if _cl[idx] == _pair[idx] and _cl[~idx] == -_pair[~idx]:
                            new_inf = _pair[idx].item()
                            dl = [-1 for i in range(len(cl.clause))]

                            all_xi = cl.clause + pair.clause
                            new_clause = pair.clause[pair.size:] + cl.clause[cl.size:] + [new_inf]

                            dl = [pair.decision_level[pair.clause.index(xi)] if xi in pair.clause else cl.decision_level[cl.clause.index(xi)] for xi in new_clause]
                            dl[-1] = decision_level

                            new_clause = Clause(new_clause)
                            new_clause.decision_level = dl
                            new_clause.size = 0
                            new_clause.value = 1

                            return new_inf, new_clause

        return None, None
        

    def bcp(self, literal, decision_level, graph):
        conflict_clause = None
        for clause in self.formula:
            if clause.value == -1:
                self.value = -1
                conflict_clause = clause
                break

            elif clause.value == 0:
                # clause.print_info()
                assert clause.size > 0
                # Implication graph is used when the lazy clause is visited
                if clause.bcp(literal, decision_level, graph) == -1:
                    self.value = -1
                    conflict_clause = clause
                    break

            elif clause.value == 1:
                continue
                
        self.value = self.get_value()
            
        return self.value, conflict_clause
            
    def unit_propagate(self, decision_level, graph=None):
        nb_clauses = len(self.formula)
        i = 0
        conflict_clause = None
        while i< nb_clauses and self.value == 0:
            clause = self.formula[i]
            if clause.is_unit():  

                unit_literal = clause.clause[0]
                if graph is not None and (unit_literal not in graph.assigned_vars):
                    graph.add_node(unit_literal, clause, decision_level)
                is_sat, conflict_clause = self.bcp(unit_literal, decision_level, graph)

                if is_sat == -1:
                    self.value == -1
                else:
                    self.value, conflict_clause = self.unit_propagate(decision_level, graph) 
                    break
            else: 
                i += 1
        return self.value, conflict_clause

    def backtrack(self, backtrack_level, graph):
        for clause in self.formula:
            clause.restore(backtrack_level, graph)
        self.value = 0

    def repair(self, lit, graph):
        for clause in self.formula:
            clause.repair_hypothesis(lit, graph)
        self.value = 0

    def add_clause(self, clause):
        self.formula += [clause]