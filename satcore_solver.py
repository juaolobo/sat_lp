from satlp import CNFLoader
from satlp import SATasLP

import numpy as np
import sys

class SATCoreSolver():

    def __init__(self, filename=None):
        if filename:
            self.cnf_loader = CNFLoader(filename)

        self.clauses = set([tuple(c) for c in cnf_loader.clauses])
        self.solutions = []

    def verify(self, witness):

        for c in self.clauses:
            abs_c = np.abs(c) 
            idx = abs_c - 1
            sgn = c/abs_c
            res = np.max(witness[idx]*sgn)
            if res == -1:
                return False

        return True

    def search_for_core(self, from_clause, det, consumed_clauses, fixing):

        # fixing = {}
        # this is changing values that were already set sometime else
        # fixing like this without much thought right now. seems to converge better
        for i,l in enumerate(from_clause):
            if i not in fixing.keys():
                fixing[abs(l)] = det[i]

        stop = False
        while len(fixing) < n and not stop:

            neg_fixing = [k*-v for k,v in fixing.items()]

            influenced_clauses = {
                tuple(c2) for c2 in self.clauses
                    if (not {x for x in np.abs(c2)}.issubset(fixing.keys())) and
                        ({c2[0],c2[1]}.issubset(neg_fixing) or 
                            {c2[0],c2[2]}.issubset(neg_fixing) or 
                            {c2[1],c2[2]}.issubset(neg_fixing))
            }

            consumed_clauses = consumed_clauses | influenced_clauses
            for ic in influenced_clauses:

                y = {xi for xi in ic if abs(xi) in fixing.keys()}

                if len(y) == len(ic) - 1:
                    idx_set = set(ic) - set(y)
                    idx = idx_set.pop()
                    fixing[abs(idx)] = np.sign(idx).item()

                else:
                    # this happens if earlier interations already covered this clause
                    # maybe fixing one clause makes it that the satisfiability is compromised
                    # may be ideal to use a queue for fixing
                    influenced_clauses = influenced_clauses - {(ic)}

            if len(influenced_clauses) == 0:    
                stop = True

        return fixing, consumed_clauses

    def _search(self, from_clause, consumed_clauses, fixing):

        remaining_clauses = self.clauses - consumed_clauses

        if len(remaining_clauses) == 0:
            return fixing

        c = from_clause
        consumed_clauses = consumed_clauses | {c}

        sgn = np.sign(c)
        possibs = [
            [sgn[j].item() if j == i else -sgn[j].item() for j in range(len(c))] 
                for i in range(len(c))
        ]
        for p in possibs:
            
            _fixing, consumed_clauses = self.search_for_core(c, p, consumed_clauses, fixing)
            _witness = np.array([fixing[i] if i in fixing.keys() else 0 for i in range(1, self.cnf_loader.n_vars+1)])

            if self.verify(_witness):
                # merge two dicts
                # this operation overwrites fixing with the values from _fixing if duplicate keys appear
                fixing = fixing | _fixing

                if len(fixing) == self.cnf_loader.n_vars:
                    witness = np.array([fixing[i] for i in range(1,self.cnf_loader.n_vars+1)])

                    if self.verify(witness):
                        return fixing
                else:
                    c = remaining_clauses.pop()
                    return self._search(c, consumed_clauses, fixing)

    def search(self):
        for c in self.clauses:
            fixing = self._search(c, set(), {})
            if fixing:
                witness = np.array([1 if fixing[i] == 1 else 0 for i in range(1, self.cnf_loader.n_vars+1)])
                self.solutions.append(witness)
                print(witness)


def verify(clauses, witness):

    for c in clauses:
        abs_c = np.abs(c) 
        idx = abs_c - 1
        sgn = c/abs_c
        res = np.max(witness[idx]*sgn)
        if res == -1:
            return False

    return True


def search(clauses, n):
    
    for c in clauses:
        fixing = _search(c, clauses, set(), {}, n)
        if fixing:
            witness = np.array([fixing[i] for i in range(1,n+1)])
            print(fixing, verify(clauses, witness))

def _search(from_clause, clauses, consumed_clauses, fixing, n):

    remaining_clauses = clauses - consumed_clauses

    if len(remaining_clauses) == 0:
        return fixing

    c = from_clause
    consumed_clauses = consumed_clauses | {c}

    sgn = np.sign(c)
    possibs = [
        [sgn[j].item() if j == i else -sgn[j].item() for j in range(len(c))] 
            for i in range(len(c))
    ]
    for p in possibs:
        
        _fixing, consumed_clauses = search_for_core(clauses, c, p, consumed_clauses, fixing)
        _witness = np.array([fixing[i] if i in fixing.keys() else 0 for i in range(1,n+1)])

        if verify(clauses, _witness):
            # merge two dicts
            # this operation overwrites fixing with the values from _fixing if duplicate keys appear
            fixing = fixing | _fixing

            if len(fixing) == n:
                witness = np.array([fixing[i] for i in range(1,n+1)])

                if verify(clauses, witness):
                    return fixing
            else:
                c = remaining_clauses.pop()
                return _search(c, clauses, consumed_clauses, fixing, n)


def search_for_core(clauses, from_clause, det, consumed_clauses, fixing):

    # fixing = {}
    # this is changing values that were already set sometime else
    # fixing like this without much thought right now. seems to converge better
    for i,l in enumerate(from_clause):
        if i not in fixing.keys():
            fixing[abs(l)] = det[i]

    stop = False
    while len(fixing) < n and not stop:

        neg_fixing = [k*-v for k,v in fixing.items()]

        influenced_clauses = {
            tuple(c2) for c2 in clauses
                if (not {x for x in np.abs(c2)}.issubset(fixing.keys())) and
                    ({c2[0],c2[1]}.issubset(neg_fixing) or 
                        {c2[0],c2[2]}.issubset(neg_fixing) or 
                        {c2[1],c2[2]}.issubset(neg_fixing))
        }

        consumed_clauses = consumed_clauses | influenced_clauses
        for ic in influenced_clauses:

            y = {xi for xi in ic if abs(xi) in fixing.keys()}

            if len(y) == len(ic) - 1:
                idx_set = set(ic) - set(y)
                idx = idx_set.pop()
                fixing[abs(idx)] = np.sign(idx).item()

            else:
                # this happens if earlier interations already covered this clause
                # maybe fixing one clause makes it that the satisfiability is compromised
                # may be ideal to use a queue for fixing
                influenced_clauses = influenced_clauses - {(ic)}

        if len(influenced_clauses) == 0:    
            stop = True

    return fixing, consumed_clauses

filename = sys.argv[1]
# filename = "cnfs/uf20-010.cnf"

cnf_loader = CNFLoader(filename=filename)

clauses = set([tuple(c) for c in cnf_loader.clauses])
n = cnf_loader.n_vars

# search(clauses, n)
# dont fully trust this because dont know how python handles attributes during recursion
solver = SATCoreSolver(filename).search()