# A collection of Efficient Global Optimization (EGO) algorithms
## Table of Contents
* [Requirements](#Requirements)
* [Standard EGO](#Standard-Efficient-Global-Optimization-algorithm)
* [Parallel EGO](#Requirements)
* [Constrained EGO](#Requirements)
* [Multiobjecitve EGO](#Requirements)
* [References](#References)
* 
## Requirements
MATLAB 2016b and above.

## Standard Efficient Global Optimization algorithm
- **EGO_EI.m** is the main file of the standard EGO algorithm. The standard EGO algorithm is implemented according to the paper of Jones et al. (1998)[[1]](#[1]). 
- For the Kriging modeling, the Gaussian correlation function is used as the corrlation function and the constant mean is used as the trend function. **Kriging_Train.m** is the function for training the Kriging model and **Kriging_Predictor.m** is the function to make predictions. I refered some codes in the book *Engineering design via surrogate modelling: a practical guide*[^2] for the Kriging model. The MATLAB **fmincon** function is used for maximizing the likehihood function to get the estimated hyperparameters when training the Kriging model.
- The expected improvement function is maximized by a real-coded genetic algorithm[^3].



## Parallel Efficient Global Optimization algorithms
Currently, three parallel EGO algorithms have been implemented.
1. The peseudo expected improvement.
2. The multipoint expected improvement.
3. The fast multipoint expected improvement.

## References
 [1] D. R. Jones, M. Schonlau, and W. J. Welch. Efficient global optimization of expensive black-box functions. Journal of Global Optimization, 1998. 13(4): 455-492.
 A. Forrester and A. Keane. Engineering design via surrogate modelling: a practical guide. 2008, John Wiley & Sons.
 K. Deb. An efficient constraint handling method for genetic algorithms. Computer Methods in Applied Mechanics and Engineering, 2000. 186(2): 311-338.
