import numpy as np
import random
import time
import argparse, sys

from satlp import HybridSolver, BooleanSolver, SATasLPFeasibility, SATasLPOptimization, SATasLPOptimizationDual, CNFLoader

def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-m", 
        "--method",
        required=False,
        default="highs-ds",
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
    lp_solver =  SATasLPOptimizationDual

    start = time.time()
    hyb_solver = HybridSolver(filename, lp_solver, method=method, track_history=True)
    witness = hyb_solver.symmetric_opt()
    print(witness)
    stop = time.time()
    print(f"Elapsed time: {stop - start}s")
    hyb_solver.verify(witness)

    sat_solver = BooleanSolver(filename, verbose=0)
    sat_solver.solve()


    # print(witness[:hyb_solver.lp_solver.n_vars()])
    # print(witness[hyb_solver.lp_solver.n_vars():2*hyb_solver.lp_solver.n_vars()])
    # print(witness[2*hyb_solver.lp_solver.n_vars():])

