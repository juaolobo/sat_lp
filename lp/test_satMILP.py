from simpleSAT_MILP import SATasLP
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
        type=str
    )
    parser.add_argument(
        "-n", 
        "--n_vars", 
        type=int
    )
    parser.add_argument(
        "-p", 
        "--run_in_parallel", 
        action='store_true', 
        help="Use this flag to run the program in multiple threads."
    )
    parser.add_argument(
        "-t", 
        "--n_threads", 
        type=int, 
        default=None, 
        help="If --run_in_parallel is active and n_threads=None, the program will use all available logical CPUs. Default is set to 'None'."
    )

    return parser

def _worker(cmb):
    
    lp_obj = SATasLP(relaxed_vars=cmb)
    lp_obj.create_lp(filename)
    res, witness = lp_obj.solve()
    ok = lp_obj.verify(witness)
    n_relaxed = len(cmb)
    row = [res, ok, witness, cmb, n_relaxed]

    return row

if __name__ == "__main__":


    cmd_parser = create_parser()
    if len(sys.argv) < 2:
        cmd_parser.print_help()
        exit(0)

    args = cmd_parser.parse_args()

    filename = f"cnfs/{args.file}"
    no_ext = args.file[:-4]
    n_threads = args.n_threads if args.run_in_parallel else 1
    n_vars = args.n_vars

    with open(f"experiments/{no_ext}-experiments.csv", "w") as f:
        writer = csv.writer(f, delimiter=',', quoting=csv.QUOTE_ALL)

        columns = ['result','is_solution','witness','relaxed_vars', 'n_relaxed']
        writer.writerow(columns)

        all_vars = range(1, n_vars+1)
        for i in tqdm(all_vars):
            cmbs = [c for c in combinations(all_vars, i)]
            print(f"Testing combinations for (20 {i})")

            with mp.Pool(n_threads) as p:
                chunksize = round(len(cmbs)/n_threads) if n_threads else round(len(cmbs)/os.cpu_count())
                csv_list = p.map(_worker, cmbs, chunksize=chunksize)

                for csv_values in csv_list:
                    writer.writerow(csv_values)