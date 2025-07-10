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

def _worker(args):

    file, method = args

    hyb_solver_feas = HybridSolver(file, SATasLPFeasibility, method=method)
    hyb_start_feas = time()
    hyb_witness_feas = hyb_solver_feas.solve()
    hyb_stop_feas = time()

    hyb_solver_opt = HybridSolver(file, SATasLPOptimization, method=method)
    hyb_start_opt = time()
    hyb_witness_opt = hyb_solver_opt.solve()
    hyb_stop_opt = time()

    sat_solver = BooleanSolver(file, verbose=0)
    bool_start = time()
    bool_witness = sat_solver.solve()
    bool_stop = time()

    ok_hyb_feas = hyb_solver_feas.verify(hyb_witness_feas)
    ok_hyb_opt = hyb_solver_opt.verify(hyb_witness_opt)
    bool_witness = sat_solver.witness_to_linear(bool_witness)
    ok_bool = hyb_solver_feas.verify(bool_witness)

    hyb_witness_feas = [i+1 if xi == 1 else -i-1 if xi == 0 else -1 for i, xi in enumerate(hyb_witness_feas)]
    hyb_witness_opt = [i+1 if xi == 1 else -i-1 if xi == 0 else -1 for i, xi in enumerate(hyb_witness_opt)]
    bool_witness = [i+1 if xi == 1 else -i-1 if xi == 0 else -1 for i, xi in enumerate(bool_witness)]

    elapsed_hyb_feas = hyb_stop_feas - hyb_start_feas
    elapsed_hyb_opt = hyb_stop_opt - hyb_start_opt
    elapsed_bool = bool_stop - bool_start

    n_learned_hyb_feas = hyb_solver_feas.cnf_handler.learnt_clauses
    n_learned_hyb_opt = hyb_solver_opt.cnf_handler.learnt_clauses
    n_learned_bool = sat_solver.nb_learnt_clause

    linear_it_feas = hyb_solver_feas.linear_it
    linear_it_opt = hyb_solver_opt.linear_it
    boolean_it_feas = hyb_solver_feas.boolean_it
    boolean_it_opt = hyb_solver_opt.boolean_it

    name = file.split("/")[-1]
    print(f"Elapsed time for our method: {elapsed_hyb_feas} (FEAS), {elapsed_hyb_opt} (OPT); elapsed time for CDCL: {elapsed_bool}")
    print(f"Finished file {name}")
    print("-------------------------------------------------------")


    row = [
        name, 
        elapsed_hyb_feas, 
        elapsed_hyb_opt, 
        elapsed_bool,
        n_learned_hyb_feas,
        n_learned_hyb_opt,
        n_learned_bool, 
        hyb_witness_feas,
        hyb_witness_opt,
        bool_witness,
        linear_it_feas,
        linear_it_opt,
        boolean_it_feas,
        boolean_it_opt
    ]

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

    n_vars = args.n_vars

    files_dir = args.dir if args.dir[-1] != '/' else args.dir[:-1]
    method = args.method
    
    files = [f"{files_dir}/{f}" for f in os.listdir(files_dir)]

    worker_fn = _worker
    experiments_file = f"experiments/data/time/uf{n_vars}-{method}.csv"

    with open(experiments_file, "w") as f:
        writer = csv.writer(f, delimiter=',', quoting=csv.QUOTE_ALL)

        columns = [
            'name', 
            'elapsed_hyb_feas', 
            'elapsed_hyb_opt', 
            'elapsed_bool',
            'n_learned_hyb_feas',
            'n_learned_hyb_opt',
            'n_learned_bool', 
            'hyb_witness_feas',
            'hyb_witness_opt', 
            'bool_witness',
            'linear_it_feas',
            'linear_it_opt',
            'boolean_it_feas',
            'boolean_it_opt'
        ]
        writer.writerow(columns)

        with mp.Pool(n_processes) as p:
            chunksize = round(len(files)/n_processes)
            chunksize = chunksize if chunksize > 0 else len(files)
            csv_list = p.map(worker_fn, zip(files, [method]*len(files)), chunksize=chunksize)
            for csv_values in csv_list:
                writer.writerow(csv_values)
