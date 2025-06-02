<h1>Quantum utils</h1>
<ul>
  <li>
    <h2>Blocking Analysis</h2>
    <p>When dealing with time series data from Monte Carlo simulations or physical measurements, data may be (too often!) statistically correlated: this implies that each new configuration of the system depends on the previous ones, and thus the system <em>remembers</em> its previous state. In these cases, the true statistical error is found to be greater than the simple standard deviation of the mean, given by: </p>
    <br>
    <p align="center">
      <img src="https://latex.codecogs.com/png.latex?\dpi{150} $\sigma_{mean} = \frac{\sigma}{\sqrt{N}}$" />
    </p>
    <br>
    <p>The correlation is characterized by a correlation time, $\tau$, which estimates how many <em>steps</em> are necessary for the system to <em>forget</em> its initial state. For data with exponential correlation: </p>
  \
  C(t) = $e^{-\frac{t}{\tau}}$ \
  \
  the theoretical statistical inefficiency is given by: \
  \
  s = 1 + $2\tau$ \
  \
  , where $\tau$ is the integrated correlation time. \
  Consecutive data points are grouped into blocks of size b, and for each block 
  </li>
- ## Density of states

  </ul>
