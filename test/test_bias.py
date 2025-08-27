from satlp import HybridSolver, BooleanSolver, SATasLPOptimizationDual
from tqdm import tqdm
import csv
import sys
import multiprocessing as mp 
import os
import argparse
from time import time
import numpy as np

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
        "-d", 
        "--dir",
        required=True,
        type=str,
        help="Directory for cnf files."
    )

    parser.add_argument(
        "-n", 
        "--n_vars", 
        required=True,
        type=int
    )

    parser.add_argument(
        "-np", 
        "--n_processes", 
        type=int, 
        default=None, 
        help="If --run_in_parallel is active and n_processes=None, the program will use all available logical CPUs. Default is set to 'None'."
    )

    parser.add_argument(
        "-s", 
        "--solver", 
        type=str, 
        default='bool', 
        help="-s [bool | feas | opt | weak_projection]."
    )

    parser.add_argument(
        "-nf", 
        "--n_formula", 
        type=int, 
        default=100, 
    )

    return parser

def _worker(args):

    file, method, solver_class, generate_cut = args

    if generate_cut is not None:
        solver = solver_class(file, lp_solver=SATasLPOptimizationDual, method=method)
        if generate_cut == "feas":
            witness = solver.optimize(generate_cut=solver.generate_feas_cut, track_history=True)

        elif generate_cut == "opt":
            witness = solver.optimize(generate_cut=generate_cut_symm, track_history=True)

        elif generate_cut == "weak_projection":
            witness = solver.optimize(generate_cut=generate_cut_via_weak_projection, track_history=True)


    else:
        solver = solver_class(file)
        witness = solver.solve()

    name = file.split("/")[-1]

    if generate_cut is not None:
        cuts = [solver.cut_to_witness(cut) for cut in solver.history]

    else:
        cuts = solver.history


    return [(name, cut) for cut in cuts]


if __name__ == "__main__":


    cmd_parser = create_parser()
    if len(sys.argv) < 2:
        cmd_parser.print_help()
        exit(0)

    args = cmd_parser.parse_args()

    n_processes = args.n_processes
    if not n_processes:
        n_processes = os.cpu_count()

    n_vars = args.n_vars

    files_dir = args.dir if args.dir[-1] != '/' else args.dir[:-1]
    method = args.method
    solver = args.solver

    if solver != "bool":    
        solver_class = HybridSolver
        if solver == "feas":
            generate_cut = solver_class.generate_feas_cut

        elif solver == "opt":
            generate_cut = solver_class.generate_cut_symm

        elif solver == "weak_projection":
            generate_cut = solver_class.generate_cut_via_weak_projection
    
    else:
        solver_class = BooleanSolver
        generate_cut = None

    n_formulas = args.n_formulas
    
    files = [f"{files_dir}/{f}" for f in os.listdir(files_dir)]
    files = np.random.choice(files, n_formulas)

    worker_fn = _worker
    experiments_file = f"experiments/data/time/uf{n_vars}-{method}.csv"

    with open(experiments_file, "w") as f:
        writer = csv.writer(f, delimiter=',', quoting=csv.QUOTE_ALL)

        columns = [
            'name', 
        ]
        writer.writerow(columns)

        with mp.Pool(n_processes) as p:
            chunksize = round(len(files)/n_processes)
            chunksize = chunksize if chunksize > 0 else len(files)
            csv_list = p.map(
                worker_fn, 
                zip(
                    files, 
                    [method]*len(files), 
                    [solver_class]*len(files),
                    [solver]*len(files)
                ),
                chunksize=chunksize
            )
            print(csv_list)
            for csv_values in csv_list:
                writer.writerow(csv_values)
