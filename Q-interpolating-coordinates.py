"""
    The present script calculates the intermediate coordinates for each atom i
    according to the interpolation formula:

        Q_i(λ) = (1 - λ) Q_i(1) + λ Q_i(2),

    and generates a new .xyz file containing the interpolated structure.

    This approach is particularly useful for studies involving polaron hopping
    or structural interpolation between two configurations.
"""

import argparse
import re
import numpy as np

ATOM_PATTERN = re.compile(r"([A-Za-z]+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)")

def calculate_shift(q1: np.ndarray, q2: np.ndarray, l: float) -> np.ndarray:
    ''' Returns q3, the interpolated position among q1 and q2. '''
    assert 0 < l < 1
    return l*q2 + (1-l)*q1

def retrieve_atom_number_per_file(filename: str) -> int:
    ''' Reads the first integer number in the given file (eg, number of total atoms). '''
    with open(filename, "r") as ifile:
        for line in ifile:
            if number := re.search(r"(\d+)", line.strip()):
                return int(number.group(1))

def interpolate_coordinates(filename1: str, filename2: str, l: float) -> None:
    ''' Interpolates the atom coordinates among two given structures (see top of the script for the formula). '''
    assert 0 < l < 1
    atom_number_file1 = retrieve_atom_number_per_file(filename1)
    atom_number_file2 = retrieve_atom_number_per_file(filename2)
    assert atom_number_file1 == atom_number_file2
    outfile = f"pol-hop-{str(l).replace(".", "-")}.xyz"
    with open(filename1, "r") as ifile1, open(filename2, "r") as ifile2, open(outfile, "w") as ofile:
        ofile.write(f"{atom_number_file1}\n\n")
        for line1, line2 in zip(ifile1, ifile2):
            if (found1 := re.search(ATOM_PATTERN, line1)) and (found2 := re.search(ATOM_PATTERN, line2)):
                ''' Checking if atoms are the same... '''
                if found1.group(1) == found2.group(1):
                    atom_symbol = found1.group(1)
                    q1 = np.array([float(found1.group(2)), float(found1.group(3)), float(found1.group(4))])
                    q2 = np.array([float(found2.group(2)), float(found2.group(3)), float(found2.group(4))])
                    q3 = calculate_shift(q1, q2, l)
                    ofile.write(f"{atom_symbol}\t{q3[0]:.8f}  {q3[1]:.8f}  {q3[2]:.8f}\n")
                else:
                    raise BufferError(f"Atoms not aligned in {filename1} and {filename1}.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str, help="The input filename containing the cartesian coordinates Q1.")
    parser.add_argument("-f", type=str, help="The input filename containing the cartesian coordinates Q2.")
    parser.add_argument("-l", type=float, help="The parameter λ for calculating Q_i(λ).")
    args = parser.parse_args()
    try:
        interpolate_coordinates(args.i, args.f, args.l)
    except FileNotFoundError as e:
        print(f"File not found: {e}")
    except Exception as e:
        print(f"An error occurred: {e}")