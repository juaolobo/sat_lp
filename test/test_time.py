from satlp import HybridSolver, BooleanSolver, SATasLPFeasibility, SATasLPOptimization
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

    hyb_solver = HybridSolver(file, SATasLPFeasibility, method=method)
    hyb_start = time()
    hyb_witness = hyb_solver.solve()
    hyb_stop = time()
    sat_solver = BooleanSolver(file, verbose=0)
    bool_start = time()
    bool_witness = sat_solver.solve()
    bool_stop = time()

    ok_hyb = hyb_solver.verify(hyb_witness)
    bool_witness = sat_solver.witness_to_linear(bool_witness)
    ok_bool = sat_solver.verify(bool_witness)

    hyb_witness = [i+1 if xi == 1 else -i-1 if xi == 0 else -1 for i, xi in enumerate(hyb_witness)]
    bool_witness = [i+1 if xi == 1 else -i-1 if xi == 0 else -1 for i, xi in enumerate(bool_witness)]

    elapsed_hyb = hyb_stop - hyb_start
    elapsed_bool = bool_stop - bool_start
    name = file.split("/")[-1]
    print(f"Elapsed time for our method: {elapsed_hyb}; elapsed time for CDCL: {elapsed_bool}")
    print(f"Finished file {name}")
    print("-------------------------------------------------------")
    row = [name, elapsed_hyb, elapsed_bool, n_learned_hyb, n_learned_bool, hyb_witness, bool_witness]

    return row


def _worker_opt(args):

    file, method = args

    hyb_solver = HybridSolver(file, SATasLPOptimization, method=method)
    hyb_start = time()
    hyb_witness = hyb_solver.solve()[:hyb_solver.lp_solver.n_vars()]
    hyb_stop = time()
    sat_solver = BooleanSolver(file, verbose=0)
    bool_start = time()
    bool_witness = sat_solver.solve()
    bool_stop = time()

    ok_hyb = hyb_solver.verify(hyb_witness[:hyb_solver.lp_solver.n_vars()])
    bool_witness = sat_solver.witness_to_linear(bool_witness)
    ok_bool = hyb_solver.verify(bool_witness)

    hyb_witness = [i+1 if xi == 1 else -i-1 if xi == 0 else -1 for i, xi in enumerate(hyb_witness)]
    bool_witness = [i+1 if xi == 1 else -i-1 if xi == 0 else -1 for i, xi in enumerate(bool_witness)]

    n_learned_bool = sat_solver.nb_learnt_clause 
    n_learned_hyb = hyb_solver.cnf_handler.learnt_clauses
    elapsed_hyb = hyb_stop - hyb_start
    elapsed_bool = bool_stop - bool_start
    name = file.split("/")[-1]

    print(f"Elapsed time for our method: {elapsed_hyb}; elapsed time for CDCL: {elapsed_bool}")
    print(f"Finished file {name}")
    print("-------------------------------------------------------")

    row = [name, elapsed_hyb, elapsed_bool, n_learned_hyb, n_learned_bool, hyb_witness, bool_witness]

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
        experiments_file = f"experiments/data/time/uf{n_vars}-feasibility-{method}.csv"

    elif lp_type == "optimization":
        worker_fn = _worker_opt
        experiments_file = f"experiments/data/time/uf{n_vars}-optimization-{method}.csv"

    with open(experiments_file, "w") as f:
        writer = csv.writer(f, delimiter=',', quoting=csv.QUOTE_ALL)

        columns = ['name', 'elapsed_hyb','elapsed_bool', 'n_learned_hyb', 'n_learned_bool', 'witness_hyb', 'witness_hyb']
        writer.writerow(columns)

        with mp.Pool(n_processes) as p:
            chunksize = round(len(files)/n_processes)
            chunksize = chunksize if chunksize > 0 else len(files)
            csv_list = p.map(worker_fn, zip(files, [method]*len(files)), chunksize=chunksize)
            for csv_values in csv_list:
                writer.writerow(csv_values)
