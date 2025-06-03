#!/usr/bin/env python

import numpy as np
from utils import create_parser
# import dpll_solver
from cdcl_solver import CDCL_Solver
import sys
import argparse

def main():
    cmd_parser = create_parser()
    if len(sys.argv) < 2:
        cmd_parser.print_help()
        exit(0)

    args = cmd_parser.parse_args()

    input_cnf_file = args.input
    verbose = args.verbose

    solver = CDCL_Solver(input_cnf_file, verbose)
    solver.solve()

if __name__ == '__main__':
    main()
