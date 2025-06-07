class Clause:
    
    def __init__(self, list_literal):
        assert len(list_literal) > 0
        self.clause = list_literal
        self.decision_level = [-1 for _ in self.clause]
        self.value = 0 # 0 = UNASSIGNED, 1 =  TRUE, -1 = FALSE
        self.size = len(self.clause)

    def print_info(self):
        print('[C] Remaining clause: ', self.clause[:self.size])
        print('[C] Truth value: ', self.value)
        print('[C] Full clause ', self.clause)
        print('[C] Details on decision_level: ', self.decision_level)

    def is_unit(self):
        if self.size == 1:
            return 1
        else:
            return 0

    def update(self, graph):
        assigned_vars = graph.assigned_vars
        sat_dl = []

        for i in range(len(self.clause)):
            # If there is one true literal in list of assigned variables => SAT and break
            if self.clause[i] in assigned_vars:

                self.decision_level[i] = graph.graph[self.clause[i]][1]
                sat_dl.append(graph.graph[self.clause[i]][1])
            # Else if the literal is false, update info
            elif -self.clause[i] in assigned_vars and self.decision_level[i] == -1:
                # if -1 in self.clause and 17 in self.clause and 9 in self.clause:
                #     breakpoint()

                self.decision_level[i] = graph.graph[-self.clause[i]][1]
            # Else continue checking next literal
            else: continue 
        
        if len(sat_dl) > 0:
            self.value = 1
            for i in range(len(self.decision_level)):
                if self.decision_level[i] == -1:
                    self.decision_level[i] = min(sat_dl) 
                elif self.decision_level[i] > min(sat_dl):
                    self.decision_level[i] = min(sat_dl)

        self.size = self.decision_level.count(-1)
        if self.size == 0 and len(sat_dl) == 0:
            self.value = -1
        

        self.clause = [x for _,x in sorted(zip(self.decision_level,self.clause), reverse=True)]
        self.decision_level.sort(reverse=True)
        self.clause = self.clause[-self.size:] + self.clause[:-self.size]
        self.decision_level = self.decision_level[-self.size:] + self.decision_level[:-self.size]


    def check_update(self, graph):

        assigned_vars = graph.assigned_vars
        sat_dl = []
        for i in range(len(self.clause)):
            # If there is one true literal in list of assigned variables => SAT and break
            if self.clause[i] in assigned_vars:
                self.decision_level[i] = graph.graph[self.clause[i]][1]
                sat_dl.append(graph.graph[self.clause[i]][1])

            # Else if the literal is false, update info
            elif -self.clause[i] in assigned_vars and self.decision_level[i] == -1:
                self.decision_level[i] = graph.graph[-self.clause[i]][1]

            # Else continue checking next literal
            else: continue 
        
        if len(sat_dl)>0:
            for i in range(len(self.decision_level)):
                if self.decision_level[i] == -1:
                    self.decision_level[i] = min(sat_dl) 
                elif self.decision_level[i] > min(sat_dl):
                    self.decision_level[i] = min(sat_dl)
            self.value = 1
            self.size = 0 

        self.size = self.decision_level.count(-1)
        ## Arrange the decision level and according literals
        self.clause = [x for _,x in sorted(zip(self.decision_level,self.clause), reverse=True)]
        self.decision_level.sort(reverse=True)

        assert self.size >= 0

        if self.size > 0:
            ## But unassigned literals now are at the end (decision_level = -1) 
            ## So move them to the head
            self.clause = self.clause[-self.size:] + self.clause[:-self.size]
            self.decision_level = self.decision_level[-self.size:] + self.decision_level[:-self.size]
            assert self.value == 0
        else:
            if self.value != 1:
                self.value = -1

        assert self.size == self.decision_level.count(-1)

    def bcp(self, literal, decision_level, graph):
        assert self.size >= 0
        assert self.size == self.decision_level.count(-1)
        # Case 1: size == 0, all literals are assigned ! check its value
        if self.size == 0: 
            assert self.value != 0

        # Case 2: If clause is unit
        elif self.size == 1:
            assert self.value == 0

            if literal == self.clause[0]: #SAT
                self.decision_level[0] = decision_level
                self.size = 0
                self.value = 1

            elif -literal == self.clause[0]: #UNSAT
                self.decision_level[0] = decision_level
                self.size = 0 
                self.value = -1
            else:
                pass

        # Case 3: clause is not unit => we need to bcp
        elif self.size > 1:

            assert self.value == 0
            self.check_update(graph)

        assert self.size == self.decision_level.count(-1)

        return self.value
            

    def restore(self, level, graph):

        self.update(graph)

        offset = 0
        for lvl in self.decision_level:
            if lvl == -1 or lvl > level:
                offset += 1

        self.size = offset
        if self.size > 0:

            self.value = 0

            self.decision_level[:self.size] = [-1 for _ in range(self.size)]

        assert self.size == self.decision_level.count(-1)

        assert max(self.decision_level) <= level

    def repair_hypothesis(self, lit, graph):

        _lit = lit if lit in self.clause else -lit if -lit in self.clause else 0
        if _lit != 0:
            self.update(graph)
            idx = self.clause.index(_lit)
            self.decision_level[idx] = -1
            self.size = self.decision_level.count(-1)
            ordered = sorted(zip(self.decision_level,self.clause), reverse=False)
            self.clause = [x[1] for x in ordered]
            self.decision_level = [x[0] for x in ordered]
            self.value = 0

    def literal_at_level(self, lvl):
        res = []
        for i in range(len(self.clause)):
            if self.decision_level[i] == lvl:
                res.append(self.clause[i])
        return res

    def get_backtrack_level(self):

        m1, m2 = -1, -1
        for x in list(set(self.decision_level)):
            if x >= m1:
                m1, m2 = x, m1
            elif x > m2:
                m2 = x
        if m2 == -1:
            m2 = m1 - 1
        return m2

    def resolution_operate(self, other, literal):

        assert (literal in self.clause and -literal in other.clause) # or  (-literal in self.clause and literal in other.clause)
        index_literal = self.clause.index(literal)
        res = self.clause[:index_literal] + self.clause[index_literal+1:]
        dl = self.decision_level[:index_literal] + self.decision_level[index_literal+1:]
        for i,l in enumerate(other.clause):
            if (abs(l) != abs(literal)):
                if l not in res and -l not in res:
                    res.append(l)
                    dl.append(other.decision_level[i])
                elif l not in res and -l in res:
                    continue
                    self.print_info()
                    print("and")
                    other.print_info()

                    res = []
                    dl = []
                    break
                elif l in res:
                    continue

        if len(res) == 0:
            return None

        resolved_clause = Clause(res)
        resolved_clause.set_decision_levels(dl)
        resolved_clause.size = resolved_clause.decision_level.count(-1)
        assert literal not in resolved_clause.clause
        assert -literal not in resolved_clause.clause
        return resolved_clause

    def set_decision_levels(self, decision_level):
        self.decision_level = decision_level
        self.clause = [x for _,x in sorted(zip(self.decision_level,self.clause), reverse=True)]
        self.decision_level.sort(reverse=True)
