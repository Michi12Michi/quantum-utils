import argparse
from collections import deque
import numpy as np
import re
from time import perf_counter

PATTERN = re.compile(r"(-?\d+\.\d+E[-+]?\d+)")
def evaluate_ipr(filename: str, pattern: str = PATTERN) -> float:
    ''' Reads the wavefunction coefficients from a .cube file and evaluates the Inverse Participation Ratio. 
        It uses re.findall, because multiple matches may occur in a single line.
        re.findall returns a list, which needs to be scanned!
    '''
    coefficient_list: deque = deque()
    with open(filename, "r") as f:
        for line in f:
            found_list: list = re.findall(pattern, line)
            if found_list:
                for coeff in found_list:
                    coefficient_list.append(float(coeff))
        np_coefficient_array: np.array = np.array(coefficient_list)
        numerator: np.float64 = np.sum(np_coefficient_array**4)
        denominator: np.float64 = np.sum(np_coefficient_array**2)**2
        return np_coefficient_array.shape[0]*(numerator/denominator)

if __name__ == "__main__":
    tic = perf_counter()
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str)
    args = parser.parse_args()
    print(evaluate_ipr(f"{args.i}"))
    toc = perf_counter()
    print(f"Process completed in {toc - tic:.2f}s.")