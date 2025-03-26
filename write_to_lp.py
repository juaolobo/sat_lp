from satlp import CNFLoader
from satlp import (
    SATasLPOptimization,
    SATasLPOptimizationEpiGraph,
    SATasLPFeasibility,
    SATasMILPOptimization,
    SATasMILPFeasibility
)
import numpy as np
import argparse
import sys

def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-f", 
        "--file",
        required=True,
        type=str,
    )

    parser.add_argument(
        "-o", 
        "--output_file",
        required=True,
        type=str,
    )

    return parser

if __name__ == "__main__":

    cmd_parser = create_parser()
    if len(sys.argv) < 2:
        cmd_parser.print_help()
        exit(0)

    args = cmd_parser.parse_args()

    filename = args.file
    output_file = args.output_file
    fixing={}
    lp = SATasLPOptimization(fixing=fixing, filename=filename, method='highs-ipm')
    lp.create_lp()

    lines = ["Minimize"]
    m, n = lp.A_ub.shape

    obj = "obj: "
    for i in range(len(lp.c)):
        if lp.c[i] > 0:
            obj += f" + {lp.c[i]} x_{i+1}"
        elif lp.c[i] < 0:
            obj += f" {lp.c[i]} x_{i+1}"

    lines.append(obj)
    lines.append("Subject To")
    for i in range(m):
        clause = f"c_{i+1}: "
        for j in range(n):
            if lp.A_ub[i][j] > 0:
                clause += f" + {lp.A_ub[i][j]} x_{j+1}"
            elif lp.A_ub[i][j] < 0:
                clause += f" {lp.A_ub[i][j]} x_{j+1}"

        clause += f" <= {lp.y_ub[i]}"
        lines.append(clause)

    if isinstance(lp.A_eq, np.ndarray):
        for i in range(lp.n_vars()):
            clause = f"c_{i+1+m}: "
            for j in range(n):
                if lp.A_eq[i][j] > 0:
                    clause += f" + {lp.A_eq[i][j]} x_{j+1}"
                elif lp.A_eq[i][j] < 0:
                    clause += f" {lp.A_eq[i][j]} x_{j+1}"

            clause += f" = {lp.y_eq[i]}"
            lines.append(clause)

    lines.append("Bounds")
    for i in range(len(lp.c)):
        if i < lp.n_vars():
            bound = f"0 <= x_{i+1} <= 1"
        
        else:
            bound = f"0 <= x_{i+1} <= 0.5"

        lines.append(bound)

    lines.append("End")
    lp_text = ("\n").join(lines)
    print(lp_text)

    with open(output_file, "w") as f:
        f.write(lp_text)


