Contains many tools useful for copula modeling in Matlab that do not exist directly in the Statistics and Machine Learning toolbox.  Highlights are:
  1) Smooth empirical copula density estimation via Beta-Kernels.
  2) Sampling from calculated empirical copula
  3) Clayton/Frank/Gumbel copula PDF and sampling for D>=2
  4) Hybrid Copula Bayesian Network construction

The (relevant) directory structure is as follows:
  - algorithms/ - contains the core copula algorithms.  Brief descriptions of the files are:
    
  | File | Description |
  | --- | --- |
  | claytoncopulapdf.m | Computes the Clayton Copula's PDF for D>=2 |
  | claytoncopularnd.m | Samples from a D>=2 Clayton Copula |
  | computeEmpiricalDiscreteProb.m | Computes empirical multinomial distribution |
  | continueRv.m | Continues realizations of a discrete RV (see http://dx.doi.org/10.1016/j.jmva.2004.01.004) |
  | empcopula_val.m | Computes value of an empirical copula at a specified point in unit hypercube |
  | empcopuladensity.m | Computes empirical copula density given pseudo-observations |
  | empcopularnd.m | Generates samples from an empirical copula |
  | estMteDensity.m | KDE with trucanted exponential distribution |
  | frankcopulapdf.m | Computes the Frank Copula's PDF for D>=2 |
  | frankcopularnd.m | Samples from a D>=2 Frank Copula |
  | gumbelcopulapdf.m | Computes the Gumbel Copula's PDF for D>=2 |
  | gumbelcopularnd.m | Samples from a D>=2 Gumbel Copula |
  | hyperFunctionError.m | Computes error between two functions in the same arbitrary dimensional plane |
  | log1mexp.m | Convenience function for log(1-exp(a)) |
  | logserrnd.m | Samples from the Log-Series distribution |
  | pseudoobs.m | Computes pseudo-observations for a given (multivariate) random vector |
  | rstable1.m | Samples from the Stable Distribution |

  - clg/ - contains code for building Conditional Linear Gaussian Bayesian Networks
    
  | File | Description |
  | --- | --- |
  | clg.m | Main class definition of CLG Model |

  - hcbn/ - contains code for building Hybrid Copula Bayesian Networks (HCBN)
    
  | File | Description |
  | --- | --- |
  | hcbn.m | Main class definition of HCBN Model |

- simulations/ - contains simulation code which uses the algorithms developed
