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
    _bootstrap = [-4, -6, 7, 9, 10, 11, -12, -14, 15, -16, -17, 18, -21, -27, -28, 29, 31, -33, 35, 37, -38, -40, 41, -43, 44, -45, 46, 47]



    solver = CDCL_Solver(input_cnf_file, verbose, _bootstrap=_bootstrap)
    solver.solve()

if __name__ == '__main__':
    main()