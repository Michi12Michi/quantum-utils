# Quantum utils
- ## Blocking Analysis
  When dealing with time series data from Monte Carlo simulations or physical measurements, data may be (too often!) statistically correlated: this implies that each new configuration of the system depends on the previous ones, and thus the system _remembers_ its previous state. In these cases, the true statistical error is found to be greater than the simple standard deviation of the mean, given by: \
  \
  $\sigma_{mean} = \frac{\sigma}{\sqrt{N}}$. \
  \
  The correlation is characterized by a correlation time, $\tau$, which estimates how many _steps_ are necessary for the system to _forget_ its initial state. For data with exponential correlation: \
  \
  C(t) = $e^{-\frac{t}{\tau}}$ \
  \
  the theoretical statistical inefficiency is given by: \
  \
  s = 1 + $2\tau$ \
  \
  , where $\tau$ is the integrated correlation time. \
- ## Density of states

  
