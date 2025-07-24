from satlp.cnf_loader import CNFLoader

from scipy.optimize import linprog
import numpy as np
from abc import ABC, abstractmethod

class SATasLPBaseclass(ABC):

    def __init__(self, filename=None, cnf_handler=None):
        self.cnf_handler = cnf_handler if cnf_handler is not None else CNFLoader(filename)

    @abstractmethod
    def _init_objects(self):
        pass

    @abstractmethod
    def _create_optimization(self):
        pass

    @abstractmethod
    def solve(self):
        pass