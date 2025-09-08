import argparse
from typing import Union
from math import sqrt
from molmass import ELEMENTS
import numpy as np
from os.path import exists
from pathlib import Path
from re import (compile, search)
from supercell_shift import SuperCellShift

ATOM_PATTERN = compile(r"([A-Za-z]+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)")

class Q:

    def __init__(
        self, 
        file_init: Union[str, Path], 
        file_final: Union[str, Path],
        constants: list
    ) -> None:

        if not exists(Path(file_init)) or not exists(Path(file_final)):
            raise FileNotFoundError("One or both of the specified files could not be found.")
        if Path(file_init).suffix != ".xyz" or Path(file_final).suffix != ".xyz":
            raise IOError("Both files must have the .xyz extension")
        if constants:
            if not all(isinstance(const, float) for const in constants):
                raise IOError("The lattice constant list must contain floating point values only.")
            else:
                self.constants = np.array(constants)
        else:
            self.constants = np.zeros(3, dtype=float)
        self.file_init = Path(file_init)
        self.file_final = Path(file_final)
        self.alignment, self.required_COM_shift = self.verify_alignment()
        self.deltaQ = self.calculate_deltaQ()
        print(f"The mass-weighted distance between the two images is: {self.deltaQ:.6f} Å * sqrt(amu).")
    
    def calculate_deltaQ(self) -> float:

        _deltaQ = 0.0
        if not self.alignment:
            print("Aligning structures according to Periodic Boundary Conditions...")
            try:
                SuperCellShift(filename=self.file_init, lattice_constants=self.constants, translation_vector=self.required_COM_shift)
            except Exception as e:
                raise RuntimeError("An error occurred during the center of mass shift. Please retry running the script with appropriate values for the lattice constants.")
            self.file_init = self.file_init.stem + "-shifted.xyz"

        with open(self.file_init, "r") as init_file, open(self.file_final, "r") as final_file:
            for i_line, f_line in zip(init_file, final_file):
                if not i_line.strip() or not f_line.strip():
                    continue
                if (found_in_i := search(ATOM_PATTERN, i_line.strip())) and (found_in_f := search(ATOM_PATTERN, f_line.strip())):
                    if found_in_i.group(1) == found_in_f.group(1):
                        mass = ELEMENTS[found_in_i.group(1)].mass
                        x_i, y_i, z_i = float(found_in_i.group(2))*sqrt(mass), float(found_in_i.group(3))*sqrt(mass), float(found_in_i.group(4))*sqrt(mass)
                        x_f, y_f, z_f = float(found_in_f.group(2))*sqrt(mass), float(found_in_f.group(3))*sqrt(mass), float(found_in_f.group(4))*sqrt(mass)
                        _deltaQ += (x_f - x_i)**2 + (y_f - y_i)**2 + (z_f - z_i)**2
                    else:
                        raise Exception("The specified .xyz files are not aligned!")
        
        return sqrt(_deltaQ)
    
    def verify_alignment(self) -> bool:
        '''
            Evaluates the center of mass of two structures and checks their physical alignment.

            Returns
            -------
            tuple[bool, np.ndarray]
                The first returned value is True if the centers of mass are aligned (within 0.01 Å), False otherwise. 
                The second returned value is a NumPy array containing the required shift to align the COMs.
        '''

        _COM1, _COM2 = np.zeros(3, dtype=float), np.zeros(3, dtype=float)
        _total_mass = 0.0
        with open(self.file_init, "r") as init_file, open(self.file_final, "r") as final_file:
            for i_line, f_line in zip(init_file, final_file):
                if not i_line.strip() or not f_line.strip():
                    continue
                if (found_in_i := search(ATOM_PATTERN, i_line.strip())) and (found_in_f := search(ATOM_PATTERN, f_line.strip())):
                    if found_in_i.group(1) == found_in_f.group(1):
                        _COM1 += ELEMENTS[found_in_i.group(1)].mass * np.array([float(found_in_i.group(2)), float(found_in_i.group(3)), float(found_in_i.group(4))])
                        _COM2 += ELEMENTS[found_in_f.group(1)].mass * np.array([float(found_in_f.group(2)), float(found_in_f.group(3)), float(found_in_f.group(4))])
                        _total_mass += ELEMENTS[found_in_f.group(1)].mass
                    else:
                        raise Exception("The specified .xyz files are not aligned!")
        _COM1 /= _total_mass
        _COM2 /= _total_mass
        _shift = _COM2 - _COM1
        _acceptable_distance_between_COMs = np.linalg.norm(_shift) <= 0.01
        print("The structures are aligned within 0.01 Å.") if _acceptable_distance_between_COMs else print("The structures are not aligned within 0.01 Å.")
        return _acceptable_distance_between_COMs, _shift

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="Q-distance-between-minima.py", 
        usage="%(prog)s -i initial_structure.xyz -f final_structure.xyz --const a_constant b_constant c_constant", 
        description="Evaluates the mass-weighted difference between the two structures defined in the .xyz files."
    )
    parser.add_argument("-i", type=str, help="")
    parser.add_argument("-f", type=str, help="")
    parser.add_argument("--const", type=float, nargs=3)

    args = parser.parse_args()
    Q(file_init=args.i, file_final=args.f, constants=args.const)
