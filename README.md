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
  </li>
  <li>
    <h2>Density of states</h2>
    <p>Manteinance... </p>
  </li>

  </ul>
