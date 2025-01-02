#import numpy as np
from ortools.linear_solver import pywraplp


def var_name( i:int ) -> str:
    return "x_" + str(i)
    
def var_pos_name( i:int ) -> str:
    return "p_" + str(i)

def var_neg_name( i:int ) -> str:
    return "n_" + str(i)

class Fml_lp:
    def __init__(self, n_vars: int, fml: list[list[int]]) -> None:
        self.n_vars   = n_vars
        self.fml      = fml
        self.vars     = None
        self.pos_vars = None
        self.neg_vars = None
        self.solver   = pywraplp.Solver.CreateSolver("GLOP")
        if not self.solver:
            raise Exception("Solver creation failed")
        
    def create_lp(self) :
        self.create_variables()
        self.create_constraints()
        self.create_optimization()
    
    def create_variables(self) -> None:
        self.vars =     [self.solver.NumVar(0.0, 1.0, var_name(i))     for i in range(1,self.n_vars+1)]
        self.pos_vars = [self.solver.NumVar(0.0, 0.5, var_pos_name(i)) for i in range(1,self.n_vars+1)]
        self.neg_vars = [self.solver.NumVar(0.0, 0.5, var_neg_name(i)) for i in range(1,self.n_vars+1)]
        print( len(self.vars), "vars,", len(self.pos_vars), "pos_vars,",len(self.neg_vars), "neg_vars" )

    def create_constraints(self):
        # noise = 1e-4
        lbs = [1 for clause in self.fml]
        for i,clause in enumerate(self.fml):
            for lit in clause:
                if lit < 0:
                    lbs[i] -= 1

        print("lbs:", lbs)
        # Clausal constraints
        constraints = [self.solver.Constraint(lbs[i], self.solver.infinity()) for i in range(len(lbs))]
        for i,clause in enumerate(self.fml):
            for lit in clause:
                var = abs(lit)
                sign = lit/(var)
                constraints[i].SetCoefficient(self.vars[var-1], sign)
        
        # var-pos-neg constraints
        for v in range(self.n_vars):
            constr = self.solver.Constraint(0.5, 0.5)
            constr.SetCoefficient(self.vars[v], 1)
            constr.SetCoefficient(self.pos_vars[v], 1)
            constr.SetCoefficient(self.neg_vars[v], -1)
            constraints.append(constr)
        
            constr = self.solver.Constraint(0.0, 0.5)
            constr.SetCoefficient(self.pos_vars[v], 1)
            constr.SetCoefficient(self.neg_vars[v], 1)
            constraints.append(constr)
        
        print("n constraints:", len(constraints))

    def create_optimization(self):
        objective = self.solver.Objective()
        lambda1 = 1/self.n_vars
        for i in range(self.n_vars):
            objective.SetCoefficient(self.pos_vars[i], 1)
            objective.SetCoefficient(self.neg_vars[i], 1)
            # Regularization
            # objective.SetCoefficient(self.vars[i], -lambda1)
        objective.SetMaximization()

    def solver(self):
        return self.solver

""" fml_teste = [[1,2],[1,-2],[-1,2],[-1,-2]]
fml_test2 = [[1,2,3],[1,-2,3],[-1,2,3],[-1,-2,3],[1,2,-3],[1,-2,-3],[-1,-2,-3],[-1,2,-3]]
fml2lp = Fml_lp(3, fml_test2)
fml2lp.create_lp()
fml2lp.solver.Solve()
print(f"Objective value = {fml2lp.solver.Objective().Value():0.9f}")
for i in range(fml2lp.n_vars):
    print(f"{var_name(i+1)} = {fml2lp.vars[i].solution_value():0.2f}", end=" ")
    print(f"{var_pos_name(i+1)} = {fml2lp.pos_vars[i].solution_value():0.2f}", end=" ")
    print(f"{var_neg_name(i+1)} = {fml2lp.neg_vars[i].solution_value():0.2f}")
 """


#msg = "Roll a dice"
#print(msg)

#print(np.random.randint(1,9))

