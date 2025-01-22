from sat_LP import SATasLPWithFixing
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
        "-nf", 
        "--n_fixed_vars", 
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

def _worker(cmb):
    ok = 0
    fixing = {abs(xi): 0 if xi < 0 else 1 for xi in cmb}
    lp_obj = SATasLPWithFixing(filename=filename, fixing=fixing)
    lp_obj.create_lp()
    status, res, witness = lp_obj.solve()
    if status != 2:
        ok = lp_obj.verify(witness)
        if ok:
            print(cmb, fixing, ok)

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
    solution_file = args.solution_file
    no_ext = args.file[:-4]
    n_vars = args.n_vars
    n_fixed_vars = args.n_fixed_vars

    with open(solution_file, "r") as f:
        solution = [int(xi) for xi in f.read().split()[1:-1]]

    combs = combinations(solution, n_fixed_vars)

    with mp.Pool(n_processes) as p:
        p.map(_worker, combs)
