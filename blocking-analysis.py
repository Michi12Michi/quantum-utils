''' 
    Usage: 
        results = BlockAnalysis(filename)
        results.print_results()
        results.create_graph()
'''

import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import re
from typing import Optional, Tuple, Dict
import warnings

class BlockAnalysis:
    """
    Implementation of blocking analysis for correlated data.
    
    This class provides the evaluation of statistical inefficiency with:
    - Automatic plateau detection
    - Uncertainty estimation via bootstrap
    - block size selection
    """

    def __init__(self, filename: str, min_blocks: int = 10, bootstrap_samples: int = 100) -> None:

        self.data: np.ndarray = np.array([])
        self.variance: float = 0.0
        self.mean: float = 0.0
        self.sigma_mean: float = 0.0
        self.sigma_mean_uncertainty: float = 0.0
        self.number_of_points: int = 0
        self.statistical_inefficiencies: np.ndarray = np.array([])
        self.inefficiency_errors: np.ndarray = np.array([])
        self.block_sizes: np.ndarray = np.array([])
        self.filename: str = filename
        self.min_blocks: int = min_blocks
        self.bootstrap_samples: int = bootstrap_samples
        self.pattern = re.compile(r"([-+]?\d+\.\d+)")
        
        # Results
        self.plateau_inefficiency: Optional[float] = None
        self.plateau_uncertainty: Optional[float] = None
        self.plateau_range: Optional[Tuple[int, int]] = None
        self.correlation_time: Optional[float] = None
        
        self._load_and_process_data()
        
    def _load_and_process_data(self) -> None:
        """Load data from file and perform initial processing."""
        self._get_data()
        if len(self.data) > 0:
            self._process_data()
            self._find_plateau()
            self._bootstrap_analysis()
            
    def _get_data(self) -> None:
        """Read data and handle errors."""
        data_list = []
        try:
            with open(self.filename, "r") as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line.startswith('#'):  # Skip empty lines and comments
                        continue
                    
                    matches = re.findall(self.pattern, line)
                    if matches:
                        # Take the first number found in each line
                        try:
                            data_list.append(float(matches[0]))
                        except ValueError:
                            warnings.warn(f"Could not parse number on line {line_num}: {line}")
                            
        except FileNotFoundError:
            raise FileNotFoundError(f"File {self.filename} not found")
        except Exception as e:
            raise RuntimeError(f"Error reading file {self.filename}: {e}")
            
        if not data_list:
            raise ValueError("No valid data found in file")
            
        self.data = np.array(data_list)
        self.number_of_points = len(self.data)
        self.mean = np.mean(self.data)
        self.variance = np.var(self.data, ddof=1)
        
        if self.variance == 0:
            raise ValueError("Data has zero variance")

    def _get_block_sizes(self) -> np.ndarray:
        """Generate optimal block sizes using logarithmic spacing."""
        max_block_size = max(10, self.number_of_points // self.min_blocks)
        
        if self.number_of_points >= 1000:
            # For large datasets, use log spacing
            block_sizes = np.unique(np.logspace(0, np.log10(max_block_size), 50).astype(int))
        else:
            # For smaller datasets, use linear spacing with finer resolution
            block_sizes = np.unique(np.linspace(1, max_block_size, min(50, max_block_size)).astype(int))
        
        valid_sizes = block_sizes[self.number_of_points // block_sizes >= self.min_blocks]
        
        return valid_sizes[valid_sizes > 0]

    def _process_data(self) -> None:
        """Process data with block size selection and error estimation."""
        block_sizes = self._get_block_sizes()
        
        if len(block_sizes) == 0:
            raise ValueError(f"Insufficient data for blocking analysis. Need at least {self.min_blocks * 10} points")
        
        inefficiencies = []
        errors = []
        valid_block_sizes = []
        
        for block_size in block_sizes:
            n_blocks = self.number_of_points // block_size
            
            if n_blocks < self.min_blocks:
                continue
                
            # Reshape data into blocks
            reshaped_data = self.data[:n_blocks * block_size].reshape((n_blocks, block_size))
            block_means = np.mean(reshaped_data, axis=1)
            block_variance = np.var(block_means, ddof=1)
            
            # Calculate statistical inefficiency
            s_estimation = block_size * block_variance / self.variance
            
            # Estimate uncertainty in s using the fact that block_variance follows chi-squared
            # Standard error of variance estimate
            se_block_var = block_variance * np.sqrt(2.0 / (n_blocks - 1))
            s_error = block_size * se_block_var / self.variance
            
            inefficiencies.append(s_estimation)
            errors.append(s_error)
            valid_block_sizes.append(block_size)
        
        self.block_sizes = np.array(valid_block_sizes)
        self.statistical_inefficiencies = np.array(inefficiencies)
        self.inefficiency_errors = np.array(errors)

    def _find_plateau(self, window_size: int = 5, tolerance: float = 0.15) -> None:
        """Automatically detect plateau region in statistical inefficiency."""
        if len(self.statistical_inefficiencies) < window_size:
            warnings.warn("Insufficient data points for plateau detection")
            return
        
        # Use a sliding window to find stable regions
        plateau_candidates = []
        
        for i in range(len(self.statistical_inefficiencies) - window_size + 1):
            window_values = self.statistical_inefficiencies[i:i + window_size]
            window_errors = self.inefficiency_errors[i:i + window_size]
            
            # Calculate relative standard deviation
            mean_val = np.mean(window_values)
            std_val = np.std(window_values)
            
            # Also consider measurement uncertainties
            avg_error = np.mean(window_errors)
            
            # A plateau should have low relative variation compared to uncertainties
            if mean_val > 0 and (std_val / mean_val < tolerance or std_val < 2 * avg_error):
                plateau_candidates.append({
                    'start_idx': i,
                    'end_idx': i + window_size - 1,
                    'mean': mean_val,
                    'std': std_val,
                    'block_size_range': (self.block_sizes[i], self.block_sizes[i + window_size - 1])
                })
        
        if plateau_candidates:
            # Choose the plateau with the largest block sizes (most reliable)
            best_plateau = max(plateau_candidates, key=lambda x: x['block_size_range'][1])
            
            self.plateau_inefficiency = best_plateau['mean']
            self.plateau_uncertainty = best_plateau['std']
            self.plateau_range = best_plateau['block_size_range']
            
            # Estimate correlation time (τ = s/2 for exponentially correlated data)
            self.correlation_time = self.plateau_inefficiency / 2.0

    def _bootstrap_analysis(self, n_bootstrap: int = None) -> None:
        """Perform bootstrap analysis to estimate uncertainty in final result."""
        if n_bootstrap is None:
            n_bootstrap = self.bootstrap_samples
            
        if self.plateau_inefficiency is None:
            return
        
        bootstrap_means = []
        bootstrap_errors = []
        
        for _ in range(n_bootstrap):
            # Resample data with replacement
            bootstrap_data = np.random.choice(self.data, size=len(self.data), replace=True)
            bootstrap_mean = np.mean(bootstrap_data)
            bootstrap_var = np.var(bootstrap_data, ddof=1)
            
            # Calculate error using plateau inefficiency
            bootstrap_error = np.sqrt(self.plateau_inefficiency * bootstrap_var / len(bootstrap_data))
            
            bootstrap_means.append(bootstrap_mean)
            bootstrap_errors.append(bootstrap_error)
        
        # Update uncertainties
        self.sigma_mean = np.sqrt(self.plateau_inefficiency * self.variance / self.number_of_points)
        self.sigma_mean_uncertainty = np.std(bootstrap_errors)

    def create_graph(self, save_path: str = None, show_plateau: bool = True) -> None:
        """Create visualization with plateau region highlighted."""
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
        plt.figure(figsize=(8, 6))
        sqrt_block_sizes = np.sqrt(self.block_sizes)
        plt.plot(sqrt_block_sizes, self.statistical_inefficiencies, 'o-', color='#526E48')
        plt.xlabel(r"$\sqrt{\tau}$")
        plt.ylabel("Statistical inefficiency")
        
        # Highlight plateau region if found
        if show_plateau and self.plateau_range is not None:
            plateau_mask = (self.block_sizes >= self.plateau_range[0]) & (self.block_sizes <= self.plateau_range[1])
            if np.any(plateau_mask):
                # Highlight plateau points
                plateau_sqrt_sizes = sqrt_block_sizes[plateau_mask]
                plateau_values = self.statistical_inefficiencies[plateau_mask]
                plt.plot(plateau_sqrt_sizes, plateau_values, 'o', color='red', markersize=8, )
                        # label=f'Plateau (s = {self.plateau_inefficiency:.2f})')
                
                # Add horizontal line for plateau value
                plt.axhline(y=self.plateau_inefficiency, color='red', linestyle='--', 
                           alpha=0.7, linewidth=2)
                
                # Add shaded region for uncertainty
                plt.axhspan(self.plateau_inefficiency - self.plateau_uncertainty,
                           self.plateau_inefficiency + self.plateau_uncertainty,
                           alpha=0.2, color='red')
                
        if save_path:
            plt.savefig(save_path)
        
        plt.show()

    def print_results(self) -> Dict:
        """Print analysis results."""
        results = {
            'mean': self.mean,
            'raw_std_error': np.sqrt(self.variance / self.number_of_points),
            'corrected_std_error': self.sigma_mean,
            'statistical_inefficiency': self.plateau_inefficiency,
            'correlation_time': self.correlation_time,
            'effective_sample_size': self.number_of_points / self.plateau_inefficiency if self.plateau_inefficiency else None
        }
        
        print("=" * 60)
        print("BLOCKING ANALYSIS RESULTS")
        print("=" * 60)
        print(f"Dataset: {self.filename}")
        print(f"Total data points: {self.number_of_points}")
        print(f"Sample mean: {self.mean:.6f}")
        print(f"Sample variance: {self.variance:.6f}")
        print()
        
        if self.plateau_inefficiency is not None:
            print("CORRELATION ANALYSIS:")
            print(f"Statistical inefficiency (s): {self.plateau_inefficiency:.3f} ± {self.plateau_uncertainty:.3f}")
            print(f"Correlation time (τ): {self.correlation_time:.3f}")
            print(f"Effective sample size: {results['effective_sample_size']:.1f}")
            print()
            
            print("ERROR ESTIMATES:")
            print(f"Naive standard error: {results['raw_std_error']:.6f}")
            print(f"Corrected standard error: {self.sigma_mean:.6f}")
            print(f"Error underestimation factor: {np.sqrt(self.plateau_inefficiency):.2f}")
            
            if self.sigma_mean_uncertainty is not None:
                print(f"Uncertainty in error estimate: ±{self.sigma_mean_uncertainty:.6f}")
        else:
            print("WARNING: No clear plateau detected!")
            print("Data may be uncorrelated or insufficient for reliable analysis.")
            print(f"Naive standard error: {results['raw_std_error']:.6f}")
        
        print("=" * 60)
        return results

    def get_final_result(self) -> Tuple[float, float]:
        """Return the final mean and its statistical error."""
        if self.sigma_mean is not None:
            return self.mean, self.sigma_mean
        else:
            # Fallback to naive estimate if plateau not found
            naive_error = np.sqrt(self.variance / self.number_of_points)
            return self.mean, naive_error