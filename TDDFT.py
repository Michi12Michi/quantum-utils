import argparse
from typing import Tuple, Union
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
from scipy.stats import norm
from pathlib import Path

plt.rcdefaults()
plt.rcParams.update({
    'axes.labelsize': 18,
    'axes.labelweight': 'bold',
    'lines.linewidth': 1,
    'lines.color': '#526E48',
    'font.family': 'Montserrat',
    'mathtext.default': 'regular' ,
    'font.size': 18,
    'legend.frameon': False,
    # 'xtick.major.size': 0,
    # 'xtick.labelsize': 0,
    'ytick.labelsize': 0,
    'ytick.major.size': 0
})

class TDDFT:
    """
    Class for processing and visualizing TDDFT excitation spectra
    from CP2K output files.

    The class reads excitation energies and oscillator strengths,
    applies Gaussian broadening, and generates a convolution spectrum.
    """

    BROADENING = 0.07

    def __init__(self, input_filename: Union[Path, str], broadening: float = BROADENING, output_filename: Union[Path, str] = None) -> None:
        """
        Initialize the TDDFT object.

        Parameters
        ----------
        input_filename : str or Path
            Path to the input file containing TDDFT excitation data.
        broadening : float, optional
            Standard deviation of the Gaussian broadening (eV).
            Default is 0.07 eV.
        output_filename : str or Path, optional
            If provided, the generated spectrum will be saved as an SVG file.

        Raises
        ------
        TypeError
            If input types are invalid.
        ValueError
            If broadening is negative.
        IOError
            If input file does not exist.
        """

        if not isinstance(input_filename, (Path, str)) or input_filename == "":
            raise TypeError("Please provide a valid input filename.")
        if (type(input_filename) == Path and not input_filename.exists()) or (type(input_filename) == str and not Path(input_filename).exists()):
                raise IOError("The input filename does not exist.")
        if broadening is None:
            broadening = self.BROADENING
        else:
            if not isinstance(broadening, float):
                raise TypeError("The broadening parameter should be a floating point value.")
            if broadening <= 0:
                raise ValueError("The broadening parameter should be a non negative value.")
        if output_filename and (not isinstance(output_filename, (Path, str)) or output_filename == ""):
            raise TypeError("Please provide a valid output filename.")
        
        self.inp = input_filename
        self.out = output_filename
        self.broad = broadening
        self.save_graph = True if output_filename else False
        self.x, self.y, self.interval, self.convolution = self.process_data()
        self.create_graph(self.save_graph)

    @staticmethod
    def gen_gauss(x, mean, std=BROADENING):
        """
        Generate a normalized Gaussian distribution.

        Parameters
        ----------
        x : np.ndarray
            X-axis values.
        mean : float
            Center of the Gaussian.
        std : float
            Standard deviation.

        Returns
        -------
        np.ndarray
            Gaussian probability density values.
        """

        return norm.pdf(x, mean, std)
    
    def process_data(self) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """
        Read excitation energies and intensities from file and
        compute the Gaussian-convoluted spectrum.

        Returns
        -------
        tuple of np.ndarray
            xTD (excitation energies), yTD (oscillator strengths), interval (energy grid), convol (convoluted spectrum)
        """
        data = np.loadtxt(fname=self.inp, usecols=(2,6), encoding="utf-8")
        xTD = data[:, 0]
        yTD = data[:, 1]
        interval = np.linspace(min(xTD) - 1, max(xTD) + 1, 1000)
        convol = np.zeros_like(interval)

        for i in range(len(xTD)):
            gaussian = self.gen_gauss(interval, xTD[i], self.broad)  
            scaled_gaussian = gaussian * yTD[i] / np.max(gaussian)  
            convol += scaled_gaussian

        return (xTD, yTD, interval, convol)
    
    def create_graph(self, save_graph: bool = False) -> None:
        """
        Plot the TDDFT spectrum.

        Parameters
        ----------
        save_graph : bool, optional
            If True, saves the figure as an SVG file.
        """

        fig, ax = plt.subplots(figsize=(8,6))
        plt.ylabel("intensity (a.u.)")
        plt.xlabel("excitation energy (eV)")

        sns.lineplot(x=self.interval, y=self.convolution, color="black", alpha=0.9)
        ax.bar(self.x, self.y, color="lime", alpha=0.6,  width=0.01)
        if save_graph:
             plt.savefig(f"{self.out}.svg")
        plt.show()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", type=str, help="The input filename containing the TDDFT data.")
    parser.add_argument("--b", type=float, help="The broadening of each gaussian curve.")
    parser.add_argument("--o", type=str, help="The output filename for the plot.", required=False)
    args = parser.parse_args()
    TDDFT(args.i, args.b, args.o)