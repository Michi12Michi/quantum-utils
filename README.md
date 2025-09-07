<h1>Quantum utils</h1>
<ul>
  <li>
    <h2>Blocking Analysis</h2>
    <p>When dealing with time series data from Monte Carlo simulations or physical measurements, data may be (too often!) statistically correlated: this implies that each new configuration of the system depends on the previous ones, and thus the system <em>remembers</em> its previous state. In these cases, the true statistical error is found to be greater than the simple standard deviation of the mean, given by: </p>
    <p align="center">
      <img src="https://latex.codecogs.com/png.latex?\dpi{150} $\sigma_{mean} = \frac{\sigma}{\sqrt{N}}$" />
    </p>
    <p>The correlation is characterized by a correlation time, $\tau$, which estimates how many <em>steps</em> are necessary for the system to <em>forget</em> its initial state.</p> 
    <p>For data with exponential correlation: </p>
    <p align="center">
      <img src="https://latex.codecogs.com/png.latex?\dpi{150} $C(t) = e^{-\frac{t}{\tau}}$" />
    </p>
    <p>, the theoretical statistical inefficiency is given by:</p>
    <p align="center">
      <img src="https://latex.codecogs.com/png.latex?\dpi{150} $s = 1 + 2\tau$" />
    </p>
    <p>, where $\tau$ is the integrated correlation time.</p>
    <p>Consecutive data points (with variance $\sigma^2_{total}$) are grouped into blocks of increasing size b, and for each block the mean and the variance are calculated, as well as the statistical inefficiency given by:</p>
    <p align="center">
      <img src="https://latex.codecogs.com/png.latex?\dpi{150} $s = b \cdot \frac{\sigma^2_{block}}{\sigma^2_{total}}$" />
    </p>
    <p>As block size b increases: </p>
    <ul>
      <li>for small sized blocks, i.e. $\rm b \ll \tau$, blocks are still correlated internally, and $\rm s < 1 + 2\tau$;</li>
      <li>for optimal sized blocks, i. e. $\rm b \approx \tau$, blocks become approximately independent, and $\rm s \approx 1 + 2\tau$;</li>
      <li>for large blocks, statistical noise increases, but the plateau value remains stable.</li>
    </ul>
    <p>
      The plateau region in the s vs. $\sqrt{b}$ plot occurs when block size becomes comparable to (or larger than) the correlation time, i.e. blocks of data become statistically independent. Therefore, the plateau value allows to determine the statistical inefficiency, the correlation time and the corrected error:
    </p>
    <p align="center">
      <img src="https://latex.codecogs.com/png.latex?\dpi{150} $\sigma_{corrected} = \sqrt{s} \cdot \sigma_{total}$" />
    </p>
    <p>The effective sample size can be determined by means of:</p>
    <p align="center">
      <img src="https://latex.codecogs.com/png.latex?\dpi{150} $\rm N_{effective} = \frac{N}{s}$" />
    </p>
      <br>
      <h5>Usage: </h5>
        <p>results = BlockAnalysis(filename)</p>
        <p>results.print_results()</p>
        <p>results.create_graph()</p>
        <p align="center">
          <img width="70%" align="center" src="https://github.com/Michi12Michi/quantum-utils/blob/main/assets/synthetic_test_results.png" />
        </p>
  </li>
      
  <li>
    <h2>Density of states (DOS)</h2>
    <p>The Density of states describes how many electronic states (both occupied and unoccupied) are available at each energy level for a system and it is crucial for understanding various material properties, including electrical conductivity, magnetic properties, and optical behavior.</p>
    <p>Mathematically, it is represented as a function D(E), where E is the energy, and it is derived from the eigenvalues of the system's Hamiltonian (the energy levels of the electrons).</p>
    <p>The Projected Density of States (pDOS) is a more detailed version of the DOS, which shows how the density of states is distributed across different atomic orbitals or specific atomic sites. While the DOS provides the total number of states at each energy level, the pDOS breaks this down further by projecting the states onto different components (such as atoms or atomic orbitals).</p>
    <p>For the DOS and the pDOS, D(E) is calculated as follows:</p>
    <p align="center">
      <img src="https://latex.codecogs.com/png.latex?\dpi{150} $\rm D(E) = \sum_{\rm i} w_{\rm i} \cdot \frac{1}{\sqrt{2 \pi \sigma^2}} exp\left(- \frac{\left(E-E_i \right)^2}{2 \sigma^2} \right)   $" />
    </p>
    <h5>Usage (1):</h5>
    <p>DOS = DensityOfStates(parameter_name: str, parameter_vector: np.ndarray, parameter_minimum: float, parameter_maximum: float, norm_factor: float, stepsize: Optional[float] = DOS_STEP, sigma: Optional[float] = DOS_SIGMA, weights: Optional[np.ndarray] = None, out_file: Union[str, Path] = None)</p>
    <p>DOS.create_graph(output_filename: Union[str, Path], parameter_vector_1: Union[np.ndarray, List, Deque], density_vector_1: Union[np.ndarray, List, Deque], curve_label1: str, x_label: str, y_label: Optional[str] = "DOS (a. u.)", parameter_vector_2: Optional[Union[np.ndarray, List, Deque]] = None, density_vector_2: Optional[Union[np.ndarray, List, Deque]] = None, curve_label2: Optional[str] = None, inverted_plot: Optional[bool] = False)</p>
    <h5>Usage (2): </h5>
    <p>DOS = DensityOfStates.from_cp2k_output(filename: Union[str, Path], subspace: str, out_file: Optional[Union[str, Path]] = None, **kwargs)</p>
    <h5>Usage (3): </h5>
    <p>pDOS = DensityOfStates.from_cp2k_pdos_output(filename: str, out_file: Optional[Union[str, Path]] = None, **kwargs)</p>
    <p align="center">
        <img width="70%" align="center" src="https://github.com/Michi12Michi/quantum-utils/blob/main/assets/test-pdos.png" />
    </p>
  </li>
   <li>
    <h2>Potential corrections</h2>
     <p>The core issue of DFT calculations is about relative energies. The Kohn-Sham eigenvalues and electrostatic potentials are only meaningful relative to some reference point within the calculation cell. Without proper referencing, it is not possible to: compare energies between different systems, calculate work functions or ionization potentials, determine band alignments at interfaces and thus connect theoretical results to experimental measurements</p>
     <p>In surface/interface calculations, the vacuum level serves as the universal energy reference because it has true physical meaning, representing the energy of a stationary electron at infinite distance from any material. Thus, connecting the reference potential to the vacuum level allows to draw experimental connections [most surface spectroscopy techniques (photoemission, inverse photoemission) naturally reference to vacuum] and ensures universality (it's the same regardless of the material system).</p>
     <h3>In progress...</h3>
   </li>
   <li>
    <h2>Supercell shift</h2>
     <p>The supercell_shift.py script allows to shift all the atoms constituting a (super)cell (with defined lattice constants) by a given translation vector, according to Periodic Boundary Conditions.</p>
     <h5>Usage: </h5>
     <p>python3 supercell_shift.py filename a_const b_const c_const x_transl y_transl z_transl</p>
   </li>
   <li>
     <h2>Mass-weighted distance between NEB images</h2>
     <p>The Q-distance-between-minima.py script allows to evaluate the mass-weighted distance between two NEB images, given their structures in aligned .xyz files. </p>
     <p>Starting from the cartesian coordinates of the two minima, <b>r<sup>A</sup></b> and <b>r<sup>B</sup></b>, for each coordinate component, k, belonging to atom a (mass m<sub>a</sub>, in amu), the mass-weighted coordinate is defined:
     <p align="center">
      <img src="https://latex.codecogs.com/png.latex?\dpi{150} $\rm q_k = \frac{r_k}{m_a}$" />
    </p>
     <p>The weighted distance between two structure is then given by:</p>
       <p align="center">
        <img src="https://latex.codecogs.com/png.latex?\dpi{150} $\rm \Delta q = \left\lVert \textbf{q}^B - \textbf{q}^A \right\rVert = \sqrt{\sum_{k=1}^{3N} \frac{(r^B_k - r^A_k)^2}{m_{ak}}}$" />
      </p>
     <h5>Usage:</h5>
     <p>python3 Q-distance-between-minima.py filename_initial filename_final</p>
   </li>
  </ul>
