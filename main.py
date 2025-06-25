import numpy as np
import random
import time
import argparse, sys

from satlp import HybridSolver, BooleanSolver, SATasLPFeasibility, SATasLPOptimization

def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-m", 
        "--method",
        required=False,
        default="highs-ipm",
        type=str,
        help="'highs-ds'= Use simplex to solve LP. 'highs-ipm' = Use interior poins method to solve LP."
    )

    parser.add_argument(
        "-f", 
        "--input_file",
        required=True,
        type=str,
        help="Input files."
    )

    parser.add_argument(
        "-t", 
        "--type",
        required=False,
        default="optimization",
        type=str,
        help="'optimization'= Use absolute value formulation. 'feasibility' = Use simple formulation without obj function."
    )

    return parser


if __name__ == "__main__":

    cmd_parser = create_parser()
    if len(sys.argv) < 2:
        cmd_parser.print_help()
        exit(0)

    args = cmd_parser.parse_args()

    filename = args.input_file
    method = args.method
    lp_solver =  SATasLPOptimization if args.type == "optimization" else SATasLPFeasibility
    start = time.time()
    hyb_solver = HybridSolver(filename, lp_solver, method=method)
    witness = hyb_solver.solve()
    stop = time.time()
    print(f"Elapsed time: {stop - start}s")
    hyb_solver.verify(witness)
