from satlp import SATasMILPFeasibility, SATasMILPOptimization
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
        help="'optimization'= Use absolute value formulation. 'feasibility'= Use simple formulation without obj function."
    )
    parser.add_argument(
        "-f", 
        "--file",
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
        "-p", 
        "--run_in_parallel", 
        action='store_true', 
        help="Use this flag to run the program in multiple threads."
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
    
    lp_obj = SATasMILPFeasibility(relaxed_vars=cmb)
    n_relaxed = len(cmb)
    lp_obj.create_lp(filename)
    s, res, witness = lp_obj.solve()
    if s == lp_obj.solver.INFEASIBLE:
        row = ["INFEASIBLE", False, [], cmb, n_relaxed]

    else:
        ok = lp_obj.verify(witness)
        row = [res, ok, witness, cmb, n_relaxed]

    return row

def _worker(cmb):
    
    lp_obj = SATasMILPOptimization(relaxed_vars=cmb)
    n_relaxed = len(cmb)
    lp_obj.create_lp(filename)
    s, res, witness = lp_obj.solve()
    if s == lp_obj.solver.INFEASIBLE:
        row = ["INFEASIBLE", False, [], cmb, n_relaxed]

    else:
        ok = lp_obj.verify(witness)
        row = [res, ok, witness, cmb, n_relaxed]

    return row

if __name__ == "__main__":


    cmd_parser = create_parser()
    if len(sys.argv) < 2:
        cmd_parser.print_help()
        exit(0)

    args = cmd_parser.parse_args()

    filename = args.file
    no_ext = args.file.split("/")[-1][:-4]
    n_processes = args.n_processes if args.run_in_parallel else 1
    if not n_processes:
        n_processes = os.cpu_count()
    n_vars = args.n_vars

    if lp_type == "feasibility":
        worker_fn = _worker_simple
        experiments_file = f"experiments/{no_ext}-milp-feasibility.csv"

    elif lp_type == "optimization":
        worker_fn = _worker_simple
        experiments_file = f"experiments/{no_ext}-milp-optimization.csv"

    with open(experiments_file, "w") as f:
        writer = csv.writer(f, delimiter=',', quoting=csv.QUOTE_ALL)

        columns = ['result','is_solution','witness','relaxed_vars', 'n_relaxed']
        writer.writerow(columns)

        all_vars = range(1, n_vars+1)
        for i in tqdm(all_vars):
            cmbs = [c for c in combinations(all_vars, i)]
            print(f"Testing combinations for (20 {i})")
            with mp.Pool(n_processes) as p:
                chunksize = round(len(cmbs)/n_processes)
                chunksize = chunksize if chunksize > 0 else len(cmbs)

                csv_list = p.map(_worker, cmbs, chunksize=chunksize)
                for csv_values in csv_list:
                    writer.writerow(csv_values)
