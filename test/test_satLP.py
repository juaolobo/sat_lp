from satlp import SATasLPFeasibility
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
        "-m", 
        "--method",
        required=True,
        default="highs-dm",
        type=str,
        help="'highs-dm'= Use simplex to solve LP. 'highs-ipm' = Use interior poins method to solve LP."
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

def _worker_ip(cmb):
    ok = 0
    fixing = {abs(xi): 0 if xi < 0 else 1 for xi in cmb}
    lp_obj = SATasLPFeasibilityIP(filename=filename, fixing=fixing, method=method)
    lp_obj.create_lp()
    x = lp_obj.solve()
    ok = lp_obj.verify(x)

    if ok:
        print(fixing)

    return ok

def _worker(cmb):
    ok = 0
    fixing = {abs(xi): 0 if xi < 0 else 1 for xi in cmb}
    lp_obj = SATasLPFeasibility(filename=filename, fixing=fixing, method=method)
    lp_obj.create_lp()
    _, _, x = lp_obj.solve()
    ok = lp_obj.verify(x)

    if ok:
        print(fixing)

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
    global method = args.method

    with open(solution_file, "r") as f:
        solution = [int(xi) for xi in f.read().split()[1:-1]]

    
    total_combs = math.comb(n_vars, n_fixed_vars)
    for worker in [_worker_ip, _worker]:
        combs = combinations(solution, n_fixed_vars)
        combs_left = total_combs
        total = 0
        print(f"Trying total of {combs_left} combinations")
        while combs_left > 0:

            n_samples = batch_size if batch_size < combs_left else combs_left
            batch = [next(combs) for i in range(n_samples)]

            with mp.Pool(n_processes) as p:
                oks = p.map(worker, batch)

            res = sum(oks)
            print(f"Convergence percentage this batch: {res/n_samples} (batch size = {n_samples})")
            total += res
            print(f"Total convergence so far: {total/total_combs}")
            combs_left -= batch_size