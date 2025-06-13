from math import (pi, sqrt,)
from scipy import constants
from collections import deque
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import re
from time import perf_counter
from typing import Deque, List, Optional, Tuple, Union

H_TO_EV = constants.physical_constants["Hartree energy in eV"][0]

class DensityOfStates:

    DOS_STEP = 0.01
    DOS_SIGMA = 0.05
    DOS_PRE_EXP_FACTOR = 1/((sqrt(2.0*pi))*DOS_SIGMA)

    def __init__(
        self,
        parameter_name: str,
        parameter_vector: np.ndarray,
        parameter_minimum: float,
        parameter_maximum: float,
        norm_factor: float,
        stepsize: Optional[float] = DOS_STEP,
        sigma: Optional[float] = DOS_SIGMA,
        weights: Optional[np.ndarray] = None,
        out_file: Optional[Union[str, Path]] = None,
    ):
        if not isinstance(parameter_name, str):
            raise ValueError("The parameter name must be a string (for example Energy, Distance, etc).")
        if not isinstance(parameter_vector, np.ndarray):
            raise ValueError("The parameter vector must be an array.")
        if not (np.issubdtype(parameter_vector.dtype, np.float64)):
            raise ValueError("The parameter vector must be an array of np.floating of float values.")
        if parameter_vector.size == 0:
            raise ValueError("The parameter vector cannot be empty.")
        if not isinstance(parameter_maximum, (float, np.float64)) or not isinstance(parameter_minimum, (float, np.float64)):
            raise ValueError("Range parameters (maximum, minimum) must be floating point values.")
        if parameter_maximum <= parameter_minimum:
            raise ValueError("Parameter maximum must be greater than parameter minimum.")
        if not isinstance(norm_factor, (float, np.float64, int)) or not isinstance(stepsize, (float, np.float64, int)) or not isinstance(sigma, (float, np.float64, int)):
            raise ValueError("Additional parameters (normalization factor, step size and sigma) must be floating point values.")
        if norm_factor <= 0:
            raise ValueError("The normalization factor must be a positive value.")
        if stepsize <= 0:
            raise ValueError("The step size must be a positive value.")
        if sigma <= 0:
            raise ValueError("The sigma value must be a positive value.")
        if weights is not None:
            if not isinstance(weights, np.ndarray):
                raise ValueError("The weights for pDOS calculation must be an array.")
            if not (np.issubdtype(weights.dtype, np.floating) or np.issubdtype(weights.dtype, np.integer)):
                raise ValueError("The weights for pDOS calculation must be an array of floating point or integer values.")
        if out_file and not isinstance(out_file, (str, Path)):
            raise ValueError("The filename for DOS/pDOS saving must be a string or a path.")
        if out_file and len(str(out_file)):
            raise ValueError("The filename for DOS/pDOS saving must be a non empty string or a path.")
        
        tic = perf_counter()
        self.parameter_name = parameter_name
        self.parameter_vector = parameter_vector
        self.parameter_minimum = parameter_minimum
        self.parameter_maximum = parameter_maximum
        self.norm_factor = norm_factor
        self.stepsize = stepsize
        self.sigma = sigma
        self.weights = weights if weights is not None else np.ones_like(self.parameter_vector)
        self.out_file = out_file
        self.fermi_energies = np.array([])
        self.density_vector, self.interval_vector = self.create_density_vector()
        toc = perf_counter()
        print(f"Execution time: {toc - tic:.3f}s")

    @classmethod
    def from_cp2k_output(cls, filename: Union[str, Path], subspace: str, out_file: Optional[Union[str, Path]] = None, **kwargs):
        """
            Create a DensityOfStates instance from a CP2K output file.
            
            Parameters:
            -----------
            filename : str or Path
                Path to the CP2K .out file.
            subspace : str
                Spin subspace identifier (e.g., "1" or "2").
            out_file : str or Path, optional
                Output file name for saving DOS data.
            **kwargs : dict
                Additional keyword arguments (stepsize, sigma, etc.).
        """

        if not isinstance(filename, (str, Path)):
            raise ValueError("The filename of the CP2K output must be a string or a path.")
        if len(str(filename)) == 0:
            raise ValueError("The filename of the CP2K output must be a non empty string or a path.")
        if isinstance(filename, Path) and not filename.exists():
            raise IOError(f"The specified file {filename} does not exist.")
        if not isinstance(subspace, str) or not len(subspace) == 1:
            raise ValueError("The subspace value must be a single character string.")
        if out_file and not isinstance(out_file, (str, Path)):
            raise ValueError("The filename for DOS/pDOS saving must be a string or a path.")
        if out_file and len(str(out_file)) == 0:
            raise ValueError("The filename for DOS/pDOS saving must be a non empty string or a path.")
        
        kohn_shams = deque()
        fermi_energies = []
        
        # State tracking for parsing
        in_occupied_section = False
        in_unoccupied_section = False
        current_spin = None
        
        patterns = {
            'occupied_header': re.compile(rf"Eigenvalues of the occupied subspace spin\s+{subspace}"),
            'unoccupied_header': re.compile(rf"Lowest Eigenvalues of the unoccupied subspace spin\s+{subspace}"),
            'separator': re.compile(r"^\s*-+\s*$"), 
            'fermi': re.compile(r"Fermi Energy \[eV\]\s*:\s*(-?\d+\.\d+)"),
            'homo_lumo': re.compile(r"HOMO - LUMO gap"),
            'force_eval_end': re.compile(r"FORCE_EVAL"),
            'convergence': re.compile(r"Reached convergence in"),
            'numbers': re.compile(r"(-?\d+.\d+)")
        }
        
        try:
            with open(filename, "r", encoding='utf-8') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    
                    if not line:
                        continue
                    
                    # Check for section headers
                    if patterns['occupied_header'].search(line):
                        in_occupied_section = True
                        in_unoccupied_section = False
                        current_spin = subspace
                        continue
                    
                    elif patterns['unoccupied_header'].search(line):
                        in_unoccupied_section = True
                        in_occupied_section = False
                        current_spin = subspace
                        continue
                    
                    # Check for section separators (lines of dashes)
                    elif patterns['separator'].search(line):
                        continue
                    
                    # Check for end of sections
                    elif (patterns['fermi'].search(line) or 
                          patterns['homo_lumo'].search(line) or
                          patterns['force_eval_end'].search(line)):
                        in_occupied_section = False
                        in_unoccupied_section = False
                        
                        # Extract Fermi energy if present
                        fermi_match = patterns['fermi'].search(line)
                        if fermi_match:
                            fermi_energies.append(float(fermi_match.group(1)))
                        continue
                    
                    # Skip convergence messages
                    elif patterns['convergence'].search(line):
                        continue
                    
                    # Extract eigenvalues if we're in the right section
                    elif (in_occupied_section or in_unoccupied_section) and current_spin == subspace:
                        # Find all numbers in the line
                        numbers = patterns['numbers'].findall(line)
                        if numbers:
                            # Convert from Hartree to eV and add to list
                            eigenvalues = [float(num) * H_TO_EV for num in numbers]
                            kohn_shams.extend(eigenvalues)
                            
        except Exception as e:
            print(f"Error reading file {filename}: {e}")
            raise

        if not kohn_shams:
            raise ValueError(f"No Kohn-Sham eigenvalues found for spin subspace {subspace}")

        stepsize = kwargs.get('stepsize', cls.DOS_STEP)
        sigma = kwargs.get('sigma', cls.DOS_SIGMA)
        
        dos = cls(
            parameter_name="Energy (eV)",
            parameter_vector=np.array(kohn_shams),
            parameter_minimum=min(kohn_shams),
            parameter_maximum=max(kohn_shams),
            norm_factor=len(kohn_shams),
            stepsize=stepsize,
            sigma=sigma,
            out_file=out_file
        )
        dos.fermi_energies = np.array(fermi_energies)
        return dos

    @classmethod
    def from_cp2k_pdos_output(cls, filename: str, out_file: Optional[Union[str, Path]] = None, **kwargs):
        """
            Creates a DensityOfStates instance from a CP2K PDOS file.
            
            Parameters
            ----------
            filename : str or Path
                Path to the CP2K .pdos file
            out_file : str or Path, optional
                Output file name for saving the DOS data
            **kwargs : dict
                Additional keyword arguments (stepsize, sigma, etc.)
        """
        if not isinstance(filename, (str, Path)):
            raise ValueError("The filename (*.pdos) of the CP2K must be a string or a path.")
        if len(str(filename)) == 0:
            raise ValueError("The filename (*.pdos) of the CP2K must be a non empty string or a path.")
        if isinstance(filename, Path) and not filename.exists():
            raise IOError(f"The specified file {filename} does not exist.")
        if out_file and not isinstance(out_file, (str, Path)):
            raise ValueError("The filename for DOS/pDOS saving must be a string or a path.")
        if out_file and len(str(out_file)) == 0:
            raise ValueError("The filename for DOS/pDOS saving must be a non empty string or a path.")
        
        squared_sums, energies = cls.pdos_energies(filename=filename)
        
        squared_sums = np.array(squared_sums)
        energies = np.array(energies)
        parameter_minimum = float(np.min(energies))
        parameter_maximum = float(np.max(energies))
        norm_factor = np.sum(squared_sums)
        stepsize = kwargs.get('stepsize', cls.DOS_STEP)
        sigma = kwargs.get('sigma', cls.DOS_SIGMA)
        
        return cls(
            parameter_name="Energy (eV)",
            parameter_vector=energies,
            parameter_minimum=parameter_minimum,
            parameter_maximum=parameter_maximum,
            norm_factor=norm_factor,
            stepsize=stepsize,
            sigma=sigma,
            weights=squared_sums,
            out_file=out_file
        )

    def create_density_vector(self,) -> tuple[np.array, np.array]:
        '''
            Calculates the normalized DOS (according to self.norm_factor) on the range between parameter_minimum and parameter_maximum. 
            Returns a tuple object of np.ndarray containing DOS and its range. If self.out_file is specified, it also produces formatted 
            text files with raw data. 
        '''
        
        bins = self.calculate_bins(
            min_value=self.parameter_minimum, 
            max_value=self.parameter_maximum, 
            stepsize=self.stepsize
            )
        density_vector = np.zeros(bins)
        interval_vector = np.linspace(self.parameter_minimum - 1.0, self.parameter_maximum + 1.0, bins)

        # first approach lol
        # for i in range(bins):
        #     density_vector[i] = np.sum(
        #         self.weights * self.DOS_PRE_EXP_FACTOR * np.exp(-((interval_vector[i] - self.parameter_vector) ** 2) / (2.0 * self.sigma ** 2))
        #     )
        energy_diff = interval_vector[:, None] - self.parameter_vector[None, :] 
        gaussian_kernel = np.exp(- (energy_diff ** 2) / (2.0 * self.sigma ** 2)) 
        density_vector = np.sum(self.weights * self.DOS_PRE_EXP_FACTOR * gaussian_kernel, axis=1)

        density_vector /= self.norm_factor
        if self.out_file:
            out_path = Path(self.out_file)
            out_path = out_path.with_suffix(".txt")
            print(f"\tWriting normalized DOS in {out_path}.txt")
            try:
                np.savetxt(f"{out_path}", np.column_stack((interval_vector, density_vector)), delimiter="\t", header=f"{self.parameter_name}\tDensity of states", comments="",)
            except Exception as e:
                print(f"An error occurred writing the file {self.out_file}.txt. {e}")
                raise
            print("\tDone!")

        return density_vector, interval_vector

    @staticmethod
    def calculate_bins(
        min_value: float, 
        max_value: float, 
        stepsize: Optional[float] = DOS_STEP,
    ) -> int:
        ''' Calculates and returns the number of bins in a given interval (with stepsize resolution). 

        Parameters
        ----------
        min_value : float or np.floating
            lower bound of the interval.
        max_value : float or np.floating
            upper bound of the interval.
        stepsize : float or np.floating or int
            defines the resolution.
        '''
        if not isinstance(min_value, (float, np.floating, int)):
            raise ValueError("Minimum range value must be a number.")
        if not isinstance(max_value, (float, np.floating, int)):
            raise ValueError("Maximum range value must be a number.")
        if min_value > max_value:
            raise ValueError("Invalid range for the binning process.")
        if not isinstance(stepsize, (float, np.floating, int)):
            raise ValueError("Stepsize value must be a number.")
        elif not stepsize > 0.0:
            raise ValueError("Stepsize must be a positive value.")
        padding = 1.5
        return int((max_value - min_value + 2*padding)/stepsize)

    @staticmethod
    def pdos_energies(filename: Union[str, Path]) -> Tuple[np.ndarray]:
        ''' 
            Reads a pdos CP2K file and calculates the square of sums foreach row (ie foreach MO). 
            Returns a tuple object of a list of squared sums of coefficients and a deque of energies.
        
            Parameters
            ----------
            filename : str or Path
                full name of the CP2K pdos file.
        '''
        
        if not isinstance(filename, (str, Path)):
            raise ValueError("The filename of the CP2K pdos must be a string or a path.")
        if len(str(filename)) == 0:
            raise ValueError("The filename of the CP2K pdos must be a non empty string or a path.")
        if isinstance(filename, Path) and not filename.exists():
            raise IOError(f"The specified file {filename} does not exist.")

        squared_sums = deque()
        energies = deque()
        try:
            data = np.loadtxt(filename, comments="#")
            if data.size == 0:
                raise ValueError(f"{filename} is empty or contains no valid data")
            energies = data[:, 1]*H_TO_EV
            kohn_shams = data[:, 3:]
            squared_sums = np.array((np.sum(kohn_shams, axis=1)**2))
        except Exception as e:
            print(f"An error occurred reading the file {filename}. {e}")
            raise
        return squared_sums, np.array(energies)

    @staticmethod
    def create_graph(
        output_filename: Union[str, Path],
        parameter_vector_1: Union[np.ndarray, List, Deque], 
        density_vector_1: Union[np.ndarray, List, Deque], 
        curve_label1: str,
        x_label: str,
        y_label: Optional[str] = "DOS (a. u.)",
        parameter_vector_2: Optional[Union[np.ndarray, List, Deque]] = None, 
        density_vector_2: Optional[Union[np.ndarray, List, Deque]] = None,
        curve_label2: Optional[str] = None,
        inverted_plot: Optional[bool] = False,
    ) -> None:
        ''' 
            Creates a plot for a DOS (eventually for the two vectors with inverted coordinates). 

            Parameters
            ----------
            output_filename : str or Path
                filename (without extension) for the .svg image plot.
            parameter_vector_1 : np.ndarray or deque or list
                x values for the first curve.
            density_vector_1 : np.ndarray or deque or list
                y values for the first curve.
            curve_label_1 : str
                descriptive label for the first curve.
            x_label : str
                descriptive label for the x axis.
            y_label : str, optional
                descriptive label for the y axis (defaults to 'DOS (a. u.)').
            parameter_vector_2 : np.ndarray or deque or list, optional
                x values for the second curve (defaults to None).
            density_vector_2 : np.ndarray or deque or list, optional
                y values for the second curve (defaults to None).
            curve_label2 : str, optional
                descriptive label for the second curve (defaults to None).
            inverted_plot : bool, optional
                flag for plotting the second curve on negative y values (defaults to False).        
        '''
        
        if not isinstance(output_filename, (str, Path)):
            raise ValueError("The filename of for the plot must be a string or a path.")
        if len(str(output_filename)) == 0:
            raise ValueError("The filename of for the plot must be a non empty string or a path.")
        output_filename = Path(output_filename)
        output_filename = output_filename.with_suffix(".svg")
        if not isinstance(parameter_vector_1, (list, deque, np.ndarray)):
            raise ValueError("The first parameter vector must be a list, or a deque or an array.")
        if not (np.issubdtype(parameter_vector_1.dtype, (np.float64, float))):   
            raise ValueError("The first parameter vector must be a list, or a deque or an array of floating point values.")
        if not isinstance(density_vector_1, (list, deque, np.ndarray)):
            raise ValueError("The first density vector must be a list, or a deque or an array.")
        if not (np.issubdtype(density_vector_1.dtype, (np.float64, float))):   
            raise ValueError("The first density vector must be a list, or a deque or an array of floating point values.")
        if not isinstance(curve_label1, str):
            raise ValueError("The label for curve (1) must be a string")
        if not isinstance(x_label, str) or len(x_label) == 0:
            raise ValueError("The label for the x axis (curve 1) must be a non empty string")
        if not isinstance(y_label, str) or len(y_label) == 0:
            raise ValueError("The label for the y axis (curve 1) must be a non empty string")
        if parameter_vector_2 is not None and not isinstance(parameter_vector_2, (deque, np.ndarray)):
            raise ValueError("The second parameter vector must be a list, or a deque or an array.")
        if parameter_vector_2 is not None and not (np.issubdtype(parameter_vector_2.dtype, (np.float64, float))):   
            raise ValueError("The second parameter vector must be a list, or a deque or an array of floating point values.")
        if density_vector_2 is not None and not isinstance(density_vector_2, (deque, np.ndarray)):
            raise ValueError("The second density vector must be a list, or a deque or an array.")
        if density_vector_2 is not None and not (np.issubdtype(density_vector_2.dtype, (np.float64, float))):   
            raise ValueError("The second density vector must be a list, or a deque or an array of floating point values.")
        if curve_label2 is not None and not isinstance(curve_label2, str):
            raise ValueError("The label for curve (2) must be a string")
        if inverted_plot and not isinstance(inverted_plot, bool):
            raise ValueError("The inverted_plot parameter must be a bool value (True|False).")

        plt.rcdefaults()
        plt.rcParams.update({
            'axes.labelsize': 18,
            'axes.labelweight': 'bold',
            'lines.linewidth': 1,
            'lines.color': '#526E48',
            'font.family': 'Montserrat',
            'font.size': 18,
            'legend.frameon': False,
            'ytick.major.size': 0,
            'ytick.labelsize': 0
        })
        fig, ax = plt.subplots(figsize=(10, 6))
        ax.plot(parameter_vector_1, density_vector_1, label=f"{curve_label1}", color='blue', alpha=0.8)
        
        if density_vector_2 is not None and parameter_vector_2 is not None:
            if inverted_plot:
                ax.plot(parameter_vector_2, -np.array(density_vector_2), label=f"{curve_label2}", color='red', alpha=0.8)
            else:
                ax.plot(parameter_vector_2, np.array(density_vector_2), label=f"{curve_label2}", color='red', alpha=0.4)

        ax.set_xlabel(f"{x_label}")  
        ax.set_ylabel(f"{y_label}")  
        ax.legend(loc="upper left", fontsize=16)
        if parameter_vector_2 is None:
            plt.xlim(np.min(parameter_vector_1), np.max(parameter_vector_1))
        else:
            plt.xlim(np.min(np.min(parameter_vector_1, parameter_vector_2)), np.max(np.max(parameter_vector_1, parameter_vector_2)))
        plt.savefig(f"{output_filename}", dpi=600)
        plt.show()