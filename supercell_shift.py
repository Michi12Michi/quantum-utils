import argparse
import numpy as np
from os.path import exists
from pathlib import Path
import re
from typing import Union

ATOM_PATTERN = re.compile(r"([A-Za-z]+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)\s*(-?\d+\.\d+)")

class SuperCellShift:
     
    def __init__(self, filename: Union[Path, str], lattice_constants: np.ndarray, translation_vector: np.ndarray):
        if filename.suffix != ".xyz" or not exists(Path(filename)):
            raise IOError("The filename parameter must be a valid, existing .xyz file.")
        if not isinstance(lattice_constants, np.ndarray) or len(lattice_constants) != 3:
            raise TypeError("The lattice_constants parameter must be a 3D Numpy array.")
        if not all(isinstance(c, float) for c in lattice_constants):
            raise TypeError("The lattice_constants parameter must be defined with floating point values.")
        if np.allclose(lattice_constants, np.zeros(3)):
            raise ValueError("The lattice_constants parameter must be defined with non zero floating point values.")
        if not isinstance(translation_vector, np.ndarray) or len(translation_vector) != 3:
            raise TypeError("The translation_vector parameter must be a 3D Numpy array.")
        if not all(isinstance(c, float) for c in translation_vector):
            raise TypeError("The translation_vector parameter must be defined with floating point values.")
        if np.allclose(translation_vector, np.zeros(3)):
            raise ValueError("The translation_vector parameter must be defined with non zero floating point values.")

        self.filename = Path(filename)
        self.lattice_constants = lattice_constants
        self.translation_vector = translation_vector
        self.atom_symbols_list: list[str] = []
        self.atom_coordinates_list: list[np.ndarray] = []

        self.shift_lattice()
        self.write_shifted_cell()

    def shift_lattice(self,) -> None:
            
        ''' Operates a lattice shift (according to PBCs), reading a xyz file and producing a xyz file. '''
        self.atom_symbols_list = []
        with open(self.filename, "r") as infile:
            for line in infile:
                line = line.strip()
                if found_atom := re.search(ATOM_PATTERN, line):
                    current_atom_replicas = []
                    atom_coordinates = np.array([float(found_atom.group(2)), float(found_atom.group(3)), float(found_atom.group(4))])
                    translated_atom = atom_coordinates + self.translation_vector
                    if self.check_coordinates(translated_atom, self.lattice_constants):
                        if not any(np.array_equal(translated_atom, coords) for coords in self.atom_coordinates_list):
                            self.atom_symbols_list.append(found_atom.group(1))
                            self.atom_coordinates_list.append(translated_atom)
                    else:
                        current_atom_replicas.extend(self.generate_periodic_replicas(atom_coordinates, self.lattice_constants))
                        for replica in current_atom_replicas:
                            translated_replica = replica + self.translation_vector
                            if self.check_coordinates(translated_replica, self.lattice_constants):
                                if not any(np.array_equal(translated_replica, coords) for coords in self.atom_coordinates_list):
                                    self.atom_symbols_list.append(found_atom.group(1))
                                    self.atom_coordinates_list.append(translated_replica)
    
    @staticmethod
    def check_coordinates(
        actual_coordinates: np.ndarray,
        lattice_constants: np.ndarray,
        ) -> bool:
        ''' The function checks if the given coordinates are within the bounds of the cell volume. 
        
            Parameters
            ----------
            actual_coordinates : np.ndarray
                3D Numpy vector containing the coordinate of a site.
            lattice_constants : np.ndarray
                3D Numpy vector containing the lattice constants of the cell.
            Returns
            -------
            bool
                Returns True if the coordinates are within the bounds of the cell volume, and False if they are outside the cell.
        '''

        if (0 <= actual_coordinates[0] <= lattice_constants[0]) and (0 <= actual_coordinates[1] <= lattice_constants[1]) and (0 <= actual_coordinates[2] <= lattice_constants[2]):
            return True
        return False
    
    def write_shifted_cell(
        self,
    ) -> None:
        
        with open(f"{self.filename.stem}-shifted.xyz", "w") as f:
            f.write(f"{len(self.atom_symbols_list)}\n\n")
            for atom in zip(self.atom_symbols_list, self.atom_coordinates_list):
                f.write(f"{atom[0]}\t{atom[1][0]:.6f}\t{atom[1][1]:.6f}\t{atom[1][2]:.6f}\n")

    @staticmethod
    def generate_periodic_replicas(
        actual_position: np.ndarray,
        lattice_constants: np.ndarray  
        ) -> list[np.array]:
            ''' Generates all periodic replicas for a given atom. '''
            periodic_replicas = []
            for dx in [-1, 0, 1]:
                for dy in [-1, 0, 1]:
                    for dz in [-1, 0, 1]:
                        if (dx == 0) and (dy == 0) and (dz == 0):
                            continue
                        periodic_replicas.append(actual_position + np.array([dx*lattice_constants[0], dy*lattice_constants[1], dz*lattice_constants[2]]))
            return periodic_replicas


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="supercell_shift.py", 
        usage="%(prog)s FILENAME a b c x y z", 
        description="Reads a cell from an .xyz file with defined constant lattices (a, b and c) and shifts it by a defined translation vector (x, y, z)."
    )
    parser.add_argument("filename", type=Path, help="")
    parser.add_argument("lattice_constants", type=float, nargs=3, metavar=('a', 'b', 'c'))
    parser.add_argument("translation_vector", type=float, nargs=3, metavar=('x', 'y', 'z'))
    args = parser.parse_args()
    SuperCellShift(
        filename=args.filename,
        lattice_constants=np.array(args.lattice_constants),
        translation_vector=np.array(args.translation_vector)
    )