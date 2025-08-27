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
        "-nf", 
        "--n_formulas", 
        type=int, 
        default=100, 
    )


    return parser

def _worker(args):

    file = args

    hyb_solver_ds = HybridSolver(file, SATasLPOptimizationDual, method='highs-ipm')
    hyb_start_ds = time()
    hyb_witness_ds = hyb_solver_ds.optimize(generate_cut=hyb_solver_ds.generate_cut_via_weak_projection)
    hyb_stop_ds = time()

    hyb_solver_ipm = HybridSolver(file, SATasLPOptimizationDual, method='highs-ds')
    hyb_start_ipm = time()
    hyb_witness_ipm = hyb_solver_ipm.optimize(generate_cut=hyb_solver_ipm.generate_cut_via_weak_projection)
    hyb_stop_ipm = time()

    ok_hyb_ds = hyb_solver_ds.verify(hyb_witness_ds)
    ok_hyb_ipm = hyb_solver_opt.verify(hyb_witness_ipm)

    hyb_witness_ds = [i+1 if xi == 1 else -i-1 if xi == 0 else -1 for i, xi in enumerate(hyb_witness_ds)]
    hyb_witness_ipm = [i+1 if xi == 1 else -i-1 if xi == 0 else -1 for i, xi in enumerate(hyb_witness_ipm)]

    elapsed_hyb_ds = hyb_stop_ds - hyb_start_ds
    elapsed_hyb_ipm = hyb_stop_ipm - hyb_start_ipm

    n_learned_hyb_ds = hyb_solver_ds.cnf_handler.learnt_clauses
    n_learned_hyb_ipm = hyb_solver_ipm.cnf_handler.learnt_clauses

    linear_it_ds = hyb_solver_ds.linear_it
    wp_it_ds = hyb_solver_ds.wp_it
    linear_it_ipm = hyb_solver_opt.linear_it
    wp_it_ipm = hyb_solver_ipm.wp_it

    name = file.split("/")[-1]
    print(f"Elapsed time for our method: {elapsed_hyb_ds} (SIMPLEX), {elapsed_hyb_ipm} (IPM)")
    print(f"Finished file {name}")
    print("-------------------------------------------------------")


    row = [
        name, 
        elapsed_hyb_ds, 
        elapsed_hyb_ipm, 
        n_learned_hyb_ds,
        n_learned_hyb_ipm,
        hyb_witness_ds,
        hyb_witness_ipm,
        linear_it_ds,
        linear_it_ipm,
        wp_it_ds,
        wp_it_ipm,
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
                zip(files),
                chunksize=chunksize
            )
            for csv_values in csv_list:
                writer.writerow(csv_values)
