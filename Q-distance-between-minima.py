import argparse
from typing import Union
from math import sqrt
from molmass import ELEMENTS
from os.path import exists
from pathlib import Path
from re import (compile, search)

ATOM_PATTERN = compile(r"([A-Za-z]+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)")

class Q:

    def __init__(
        self, 
        file_init: Union[str, Path], 
        file_final: Union[str, Path]
    ) -> None:

        if not exists(Path(file_init)) or not exists(Path(file_final)):
            raise FileNotFoundError("One or both of the specified files could not be found.")
        if Path(file_init).suffix != ".xyz" or Path(file_final).suffix != ".xyz":
            raise IOError("Both files must have the .xyz extension")
        self.file_init = Path(file_init)
        self.file_final = Path(file_final)
        self.deltaQ = self.calculate_deltaQ()
        print(f"The mass-weighted distance between the two images is: {self.deltaQ:.6f} ang / sqrt(amu).")
    
    def calculate_deltaQ(self) -> float:
        _deltaQ = 0.0

        with open(self.file_init, "r") as init_file, open(self.file_final, "r") as final_file:
            for i_line, f_line in zip(init_file, final_file):
                if not i_line.strip() or not f_line.strip():
                    continue
                if (found_in_i := search(ATOM_PATTERN, i_line.strip())) and (found_in_f := search(ATOM_PATTERN, f_line.strip())):
                    if found_in_i.group(1) == found_in_f.group(1):
                        mass = ELEMENTS[found_in_i.group(1)].mass
                        x_i, y_i, z_i = float(found_in_i.group(2))/mass, float(found_in_i.group(3))/mass, float(found_in_i.group(4))/mass
                        x_f, y_f, z_f = float(found_in_f.group(2))/mass, float(found_in_f.group(3))/mass, float(found_in_f.group(4))/mass
                        _deltaQ += (x_f - x_i)**2 + (y_f - y_i)**2 + (z_f - z_i)**2
                    else:
                        raise Exception("The specified .xyz files are not aligned!")
        
        return sqrt(_deltaQ)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="Q-distance-between-minima.py", 
        usage="%(prog)s filename_1 filename_2", 
        description="Evaluates the mass-weighted difference between the two structures defined in the .xyz files."
    )
    parser.add_argument("filename_1", type=str, help="")
    parser.add_argument("filename_2", type=str, help="")
    args = parser.parse_args()
    Q(file_init=args.filename_1, file_final=args.filename_2)
