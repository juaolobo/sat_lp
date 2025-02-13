from satlp import SATasLPFeasibility, SATasLPOptimization
from tqdm import tqdm
from itertools import combinations
import csv
import sys
import multiprocessing as mp 
import os
import argparse

def create_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-t", 
        "--type",
        required=True,
        default="optimization",
        type=str,
        help="'optimization'= Use absolute value formulation. 'feasibility' = Use simple formulation without obj function."
    )
    parser.add_argument(
        "-f", 
        "--file",
        required=True,
        type=str
    )
    parser.add_argument(
        "-sf", 
        "--solution_file",
        required=True,
        type=str
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


    return parser

def _worker_simple(cmb):
    ok = 0
    fixing = {abs(xi): 0 if xi < 0 else 1 for xi in cmb}
    n_fixed = len(fixing)
    lp_obj = SATasLPFeasibility(filename=filename, fixing=fixing)
    lp_obj.create_lp()
    status, res, witness = lp_obj.solve()

    if s == lp_obj.solver.INFEASIBLE:
        row = ["INFEASIBLE", False, [], n_fixed, fixing]

    else:
        ok = lp_obj.verify(witness)
        row = [res, ok, witness, n_fixed, fixing]
    
    return row

def _worker(cmb):
    ok = 0
    fixing = {abs(xi): 0 if xi < 0 else 1 for xi in cmb}
    n_fixed = len(fixing)
    lp_obj = SATasLPOptimization(filename=filename, fixing=fixing)
    lp_obj.create_lp()
    status, res, witness = lp_obj.solve()

    if s == lp_obj.solver.INFEASIBLE:
        row = ["INFEASIBLE", False, [], n_fixed, fixing]

    else:
        ok = lp_obj.verify(witness)
        row = [res, ok, witness, n_fixed, fixing]
    
    return row

if __name__ == "__main__":


    cmd_parser = create_parser()
    if len(sys.argv) < 2:
        cmd_parser.print_help()
        exit(0)

    args = cmd_parser.parse_args()

    n_processes = args.n_processes
    if not n_processes:
        n_processes = os.cpu_count()

    filename = args.file
    lp_type = args.type

    if lp_type == "feasibility":
        worker_fn = _worker_simple
        experiments_file = f"experiments/data/{no_ext}-simplex-feasibility.csv"

    elif lp_type == "optimization":
        worker_fn = _worker_simple
        experiments_file = f"experiments/data/{no_ext}-simplex-optimization.csv"

    solution_file = args.solution_file
    no_ext = args.file.split("/")[-1][:-4]
    n_vars = args.n_vars

    with open(solution_file, "r") as f:
        solution = [int(xi) for xi in f.read().split()[1:-1]]

    with open(experiments_file, "w") as f:
        writer = csv.writer(f, delimiter=',', quoting=csv.QUOTE_ALL)

        columns = ['result','is_solution','witness', 'n_fixed', 'fixing']
        writer.writerow(columns)

        all_vars = range(1, n_vars+1)
        for i in tqdm(all_vars):
            cmbs = [c for c in combinations(solution, i)]
            print(f"Testing combinations for ({n_vars} {i})")
            with mp.Pool(n_processes) as p:
                chunksize = round(len(cmbs)/n_processes)
                chunksize = chunksize if chunksize > 0 else len(cmbs)

                csv_list = p.map(_worker, cmbs, chunksize=chunksize)
                for csv_values in csv_list:
                    writer.writerow(csv_values)
