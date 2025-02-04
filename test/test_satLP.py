from satlp import SATasLPFeasibilityIP, SATasLPFeasibility
from tqdm import tqdm
from itertools import combinations
import csv
import sys
import multiprocessing as mp 
import os
import argparse
import math

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
    parser.add_argument(
        "-b", 
        "--batch_size", 
        type=int, 
        default=10000, 
    )


    return parser

def _worker(cmb):
    ok = 0
    fixing = {abs(xi): 0 if xi < 0 else 1 for xi in cmb}
    lp_obj = SATasLPFeasibilityIP(filename=filename, fixing=fixing)
    lp_obj.create_lp()
    x = lp_obj.solve()
    ok = lp_obj.verify(x)
    if ok:
        print(fixing, x)

    return ok

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
    batch_size = args.batch_size

    with open(solution_file, "r") as f:
        solution = [int(xi) for xi in f.read().split()[1:-1]]

    combs = combinations(solution, n_fixed_vars)
    
    combs_left = math.comb(n_vars, n_fixed_vars)

    while combs_left > 0:
        batch = [next(combs) for i in range(batch_size)]

        with mp.Pool(n_processes) as p:
            oks = p.map(_worker, batch)

        print(sum(oks)/math.comb(n_vars, n_fixed_vars))
        combs_left -= consume