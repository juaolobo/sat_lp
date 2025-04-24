import argparse

def create_parser():
    argparser = argparse.ArgumentParser(description=__doc__)
    argparser.add_argument(
        '-i', '--input',
        metavar='I',
        default='',
        help='The DIMACS file')
    argparser.add_argument(
        '-v', '--verbose',
        default=1,    
        help='Verbose option')
    return argparser