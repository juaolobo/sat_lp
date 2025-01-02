from unittest import main, TestCase
from Fml_lp import *

fml_test1 = [[1,2],[1,-2],[-1,2],[-1,-2]]
fml_test2 = [[1,2,3],[1,-2,3],[-1,2,3],[-1,-2,3],[1,2,-3],[1,-2,-3],[-1,-2,-3],[-1,2,-3]]

def frac(x:float) -> float:
    return x - int(x)


class TestFml_lp(TestCase):
    def test_2_variables_unsat(self):
        
        nvars = 2
        fml2lp = Fml_lp(nvars, fml_test1)
        fml2lp.create_lp()
        fml2lp.solver.Solve()

        print(f"Objective value = {fml2lp.solver.Objective().Value():0.9f}")
        for i in range(fml2lp.n_vars):
            print(f"{var_name(i+1)} = {fml2lp.vars[i].solution_value():0.2f}", end=" ")
            print(f"{var_pos_name(i+1)} = {fml2lp.pos_vars[i].solution_value():0.2f}", end=" ")
            print(f"{var_neg_name(i+1)} = {fml2lp.neg_vars[i].solution_value():0.2f}")

        expected = nvars/2
        result = fml2lp.solver.Objective().Value()
        self.assertEqual(result, expected)

    def test_2_variables_sat(self):
        
        nvars = 2
        fml = fml_test1[1:]
        fml2lp = Fml_lp(nvars, fml)
        fml2lp.create_lp()
        fml2lp.solver.Solve()

        print(f"Objective value = {fml2lp.solver.Objective().Value():0.9f}")
        for i in range(fml2lp.n_vars):
            print(f"{var_name(i+1)} = {fml2lp.vars[i].solution_value():0.2f}", end=" ")
            print(f"{var_pos_name(i+1)} = {fml2lp.pos_vars[i].solution_value():0.2f}", end=" ")
            print(f"{var_neg_name(i+1)} = {fml2lp.neg_vars[i].solution_value():0.2f}")

        expected = nvars/2
        result = fml2lp.solver.Objective().Value()
        self.assertEqual(result, expected)

if __name__ == '__main__':
    main()
