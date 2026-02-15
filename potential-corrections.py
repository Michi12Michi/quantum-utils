import argparse
from collections import deque
import mmap
import numpy as np
from numbers import Real
import matplotlib.pyplot as plt
from pathlib import Path
from re import (compile, findall, match)
import scipy.constants as cnst
from scipy.signal import find_peaks
from scipy.ndimage import uniform_filter1d
from typing import Tuple

patterns = {
    "general_infos": compile(r"^\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)$"),
    "coordinates": compile(r"^\s+(\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)\s+(-?\d+\.\d+)$"),
    "potentials": compile(r"\s+(-?\d+\.\d+E[+-]?\d+)")
}

H_TO_EV = cnst.physical_constants["atomic unit of electric potential"][0]
SIMILARITY_THRESHOLD = 0.01

class PotentialCorrection:

    def __init__(self, filename_path: str, c_lattice_constant: Real, fermi_level: Real) -> None:

        if len(filename_path) == 0:
            raise ValueError("The filename path must be a valid (non-null) string.")
        if not Path(filename_path).exists():
            raise IOError("The filename does not exist.")
        if Path(filename_path).suffix != ".cube":
            raise IOError("The filename has not the .cube CP2K extension.")
        try:
            c_lattice_constant = float(c_lattice_constant)
            if c_lattice_constant <= 0.0:
                raise ValueError("The c lattice constant must be a positive floating point value.")
        except ValueError as e:
            print(f"The c lattice constant must be a floating point value: {e}.")
            raise
        if fermi_level is not None:
            if not isinstance(fermi_level, float):
                raise TypeError("The Fermi level must be a floating point value.")
            self.fermi = fermi_level
        else:
            self.fermi = None
        
        self.filename_path = Path(filename_path)
        self.c_lattice_constant = float(c_lattice_constant)
        self.dx, self.dy, self.dz = 0, 0, 0
        self.potentials: np.ndarray = self.read_file()        
        self.corrected_potentials = self.autocorrect_potentials()
        self._loc_1, self._loc_2 = self.find_maxima_indexes()
        _tmp_min, _tmp_max = np.min([self._loc_1, self._loc_2]), np.max([self._loc_1, self._loc_2])
        self._loc_1, self._loc_2 = _tmp_min, _tmp_max
        self.running_averages = self.running_average(self.potentials, window_size=np.abs(self._loc_1 - self._loc_2))
        self.uncorrected_vacuum_potential = np.max(self.potentials)
        self.mean_potential = np.mean(self.running_averages[self._loc_1:self._loc_2])
        
        self.log_data()
        self.plot_potentials()
        
    def read_file(self) -> np.ndarray:
        '''
            Reads a CP2K .cube file generated with the &V_HARTREE_CUBE keyword, and extracts 
            the electrostatic potential values evaluated over a slab.

            The .cube file is expected to contain:
            - A header (first two lines), which is ignored.
            - Line 3: number of atoms in the system.
            - Lines 4-6: grid dimensions (dx, dy, dz) for the supercell.
            - Lines 7 to (6 + number_of_atoms): atomic information (ignored).
            - Remaining lines: electrostatic potential values (in Hartree units).

            Populates the following attributes:
            - self.dx, self.dy, self.dz: grid dimensions.
            - self.potentials: deque of extracted potential values (floats).

            Returns
            -------
            np.ndarray
                A Numpy array of mean electrostatic potentials (in eV) calculated on dz slices.
        '''

        _number_of_atoms = 0
        _line_counter = -1
        _potentials: deque = deque()
        try:
            with open(self.filename_path, "rb") as f:
                mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
                start = 0
                while start < len(mm):
                    end = mm.find(b'\n', start)
                    if end == -1:  
                        line = mm[start:]
                        start = len(mm)
                    else:
                        line = mm[start:end]
                        start = end + 1  
                    line_str = line.decode('utf-8')
                    _line_counter += 1
                    if _line_counter == 0 or _line_counter == 1:
                        pass
                    else:
                        if 2 <= _line_counter <= 5:
                            if _found_pattern := match(patterns["general_infos"], line_str):
                                if _line_counter == 2:
                                    _number_of_atoms = int(_found_pattern.group(1))
                                elif _line_counter == 3:
                                    self.dx = int(_found_pattern.group(1))
                                elif _line_counter == 4:
                                    self.dy = int(_found_pattern.group(1))
                                elif _line_counter == 5:
                                    self.dz = int(_found_pattern.group(1))
                            else:
                                raise ValueError(f"Line {_line_counter}: Expected general info pattern.")
                        elif _line_counter <= _number_of_atoms + 5:
                            continue
                        else:
                            if _found_pattern := findall(patterns["potentials"], line_str):
                                for _i in _found_pattern:
                                    _potentials.append(float(_i)*H_TO_EV)
                            else:
                                raise ValueError(f"Line {_line_counter}: Invalid potential format.")
        except FileNotFoundError:
            print(f"The file provided {self.filename_path} was not found.")
            raise
        except Exception as e:
            raise IOError(f"An error occurred while reading {self.filename_path}: {e}.")
        
        if len(_potentials) != self.dx*self.dy*self.dz:
            raise ValueError(f"Unexpected number of potential values: expected {self.dx*self.dy*self.dz}, got {len(_potentials)}.")

        _potentials = np.array(_potentials).reshape((self.dx, self.dy, self.dz))
        return np.mean(_potentials, axis=(0, 1))

    @staticmethod
    def running_average(values: np.ndarray, window_size: int) -> np.ndarray:
        return uniform_filter1d(values, size=window_size, mode="nearest")
    
    def log_data(self) -> None:
        '''
            Logs physical information to the console and raw data for further plotting and/or data analysis.
        '''
        print(f"(Uncorrected) Vacuum potential: {self.uncorrected_vacuum_potential} eV")
        print(f"Mean potential: {self.mean_potential} eV")
        print(f"Ionization potential: {self.uncorrected_vacuum_potential - (self.mean_potential + self.fermi)}") if self.fermi is not None else ...
        x = np.array([self.c_lattice_constant/self.dz*i for i, _ in enumerate(self.potentials)])
        data = np.column_stack((x, self.potentials, self.corrected_potentials))
        np.savetxt(f"{self.filename_path.stem}-RAWDATA.dat", data, fmt='%.6f', header="slab length z (Å)\tMean Potential (eV)\tCorrected Potential (eV)", delimiter='\t',)

    def find_maxima_indexes(self) -> Tuple[int]:
        '''
            Leveraging the scipy.signals.find_peaks function, detects all local maxima in the mean potential trend and selects two
            adjacent peaks. Their distance is used for setting a suitable window for the running average method.

            Returns
            -------
            Tuple[int]
                Tuple of indexes for the detected potential local maxima.
        '''

        _peaks = find_peaks(self.potentials)
        # construct an ordered list (key descending) of tuples (index, potential_value)
        peak_indexes_and_values = sorted([(index, self.potentials[index]) for index in _peaks[0]], key=lambda x: x[1], reverse=True)

        if not peak_indexes_and_values:
            raise ValueError("No peaks detected in the mean potential function.")
    
        _approx_vacuum = peak_indexes_and_values[0][1]
        significant_peaks = [
            item for item in peak_indexes_and_values
            if np.abs(item[1] - _approx_vacuum) > SIMILARITY_THRESHOLD
        ]

        if len(significant_peaks) < 2:
            raise ValueError("No distinct peaks found other than the vacuum level.")

        # bad, old approach :D
        # local_1, local_2 = significant_peaks[:2]
        local_1 = significant_peaks[0]
        local_2 = [
            item for item in significant_peaks
            if np.abs(item[1] - local_1[1]) > SIMILARITY_THRESHOLD
        ][0]

        return (local_1[0], local_2[0])
        
    def autocorrect_potentials(self) -> np.ndarray:
        '''
            Shifts all the mean potentials upwards, in order to assign ca. 0 eV to the vacuum level.

            Returns
            -------
            np.ndarray
                The mean potential Numpy array shifted towards the vacuum level (ca 0 eV).
        '''
        
        return self.potentials - np.max(self.potentials)
    
    def plot_potentials(self) -> None:

        plt.rcParams.update({
            'axes.labelsize': 18,
            'axes.labelweight': 'bold',
            'lines.linewidth': 1,
            'lines.color': '#526E48',
            'font.family': 'Montserrat',
            'font.size': 18,
            'legend.frameon': False,
            'ytick.major.size': 3.5
        })
        plt.rcParams.update({'mathtext.default': 'regular'})
        plt.figure(figsize=(10,6))
        plt.ylabel("energy vs vacuum (eV)")
        plt.xlabel(rf"slab length (z, $\AA$)")        
        x = np.array([self.c_lattice_constant/self.dz*i for i, _ in enumerate(self.potentials)])
        plt.xlim(min(x), max(x))
        plt.plot(x, self.potentials)
        plt.plot(x, self.running_averages, "pink", alpha=0.6, linestyle="dashdot")
        plt.scatter(x[self._loc_1], self.potentials[self._loc_1])
        plt.scatter(x[self._loc_2], self.potentials[self._loc_2])
        plt.axhline(np.mean(self.running_averages[self._loc_1:self._loc_2]), color="green", alpha=0.6, linestyle="dashed")
        plt.savefig(f"{self.filename_path.stem}.svg")
        plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str)
    parser.add_argument("-c", type=float, help="The c lattice parameter of the slab.")
    parser.add_argument("--f", type=float, help="The calculated Fermi level for the slab.", required=False)
    args = parser.parse_args()
    PotentialCorrection(args.i, args.c, args.f)