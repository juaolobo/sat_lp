from satlp import SATasLPFeasibility, SATasLPOptimization
from tqdm import tqdm
import csv
import sys
import multiprocessing as mp 
import os
import argparse
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
        "-t", 
        "--type",
        required=False,
        default="optimization",
        type=str,
        help="'optimization'= Use absolute value formulation. 'feasibility' = Use simple formulation without obj function."
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

def _worker_feas(args):


    file, method = args

    new_fixing = {}
    fixing = {}
    while True:

        lp_obj = SATasLPFeasibility(filename=file, fixing=fixing, method=method)
        lp_obj.create_lp()
        n = lp_obj.n_vars()
        witness = lp_obj.solve()
        for i, xi in enumerate(witness[:n]):
            if xi.is_integer():
                new_fixing[i+1] = xi

        if fixing != new_fixing:
            fixing = new_fixing

        else:
            break

    
    integral_vars = len(fixing)

    ok = lp_obj.verify(witness)
    name = file.split("/")[-1]
    print(f"Program achieved {integral_vars} integral variables")
    print(f"Finished file {name}")
    print("-------------------------------------------------------")

    row = [name, ok, len(fixing), witness]

    return row


def _worker_opt(args):
    
    file, method = args

    new_fixing = {}
    fixing = {}
    while True:

        lp_obj = SATasLPOptimization(filename=file, fixing=fixing, method=method)
        lp_obj.create_lp()
        n = lp_obj.n_vars()
        witness = lp_obj.solve()
        for i, xi in enumerate(witness[:n]):
            if xi.is_integer():
                new_fixing[i+1] = xi

        if fixing != new_fixing:
            fixing = new_fixing

        else:
            break

    
    integral_vars = len(fixing)
    ok = lp_obj.verify(witness)
    name = file.split("/")[-1]

    print(f"Program achieved {integral_vars} integral variables")
    print(f"Finished file {name}")
    print("-------------------------------------------------------")

    row = [name, ok, len(fixing), witness]

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

    lp_type = args.type
    n_vars = args.n_vars

    files_dir = args.dir if args.dir[-1] != '/' else args.dir[:-1]
    method = args.method
    
    files = [f"{files_dir}/{f}" for f in os.listdir(files_dir)]

    if lp_type == "feasibility":
        worker_fn = _worker_feas
        experiments_file = f"experiments/data/integral_conv/uf{n_vars}-feasibility-{method}.csv"

    elif lp_type == "optimization":
        worker_fn = _worker_opt
        experiments_file = f"experiments/data/integral_conv/uf{n_vars}-optimization-{method}.csv"

    with open(experiments_file, "w") as f:
        writer = csv.writer(f, delimiter=',', quoting=csv.QUOTE_ALL)

        columns = ['name', 'satisfied','int_vars','witness']
        writer.writerow(columns)

        with mp.Pool(n_processes) as p:
            chunksize = round(len(files)/n_processes)
            chunksize = chunksize if chunksize > 0 else len(files)
            csv_list = p.map(worker_fn, zip(files, [method]*len(files)), chunksize=chunksize)
            for csv_values in csv_list:
                writer.writerow(csv_values)
