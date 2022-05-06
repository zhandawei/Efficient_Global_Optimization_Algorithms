# A collection of Efficient Global Optimization (EGO) algorithms
## Table of Contents
* [Requirements](#Requirements)
* [Standard EGO](#Standard-Efficient-Global-Optimization-algorithm)
* [Parallel EGO](#Parallel-Efficient-Global-Optimization-algorithms)
* [Constrained EGO](#Constrained-Efficient-Global-Optimization-algorithms)
* [Multiobjecitve EGO](#Multiobjective-Efficient-Global-Optimization-algorithms)
* [References](#References)


## Requirements
MATLAB 2016b and above.

## Standard Efficient Global Optimization algorithm
1. **The standard EGO algorithm** (*EGO_EI.m*) [^1].
2. For the Kriging modeling, the Gaussian correlation function is used as the corrlation function and the constant mean is used as the trend function. I refered some codes in the book *Engineering design via surrogate modelling: a practical guide* [^2] for the Kriging model. 
3. The MATLAB **fmincon** function is used for maximizing the likehihood function to get the estimated hyperparameters when training the Kriging model. The expected improvement function is maximized by a real-coded genetic algorithm [^3].

## Parallel Efficient Global Optimization algorithms
1. **the Kriging Believer approach** (*EGO_KB.m*) [^4].
2. **the Constant Liar approach** (*EGO_CL.m*) [^4].
3. **the peseudo expected improvement** (*EGO_Pseudo_EI.m*) [^5] .
4. **the multipoint expected improvement** (*EGO_qEI.m*) [^6].
5. **the fast multipoint expected improvement** (*EGO_FqEI.m*) [^7].

## Constrained Efficient Global Optimization algorithms
1. **the constrained expected improvemment** (*EGO_Constrained_EI.m*) [^8].
2. **the pseudo constrained expected improvement** (*EGO_Pseudo_Constrained_EI.m*) [^9].


## Multiobjective Efficient Global Optimization algorithms
1. **the ParEGO (Pareto EGO)** (*ParEGO.m*) [^10].
2. **the expected improvement matrix** (*EGO_EIM_Euclidean.m*,*EGO_EIM_Hypervolume.m*,*EGO_EIM_Maximin.m*) [^11].






## References
[^1]: D. R. Jones, M. Schonlau, and W. J. Welch. Efficient global optimization of expensive black-box functions. Journal of Global Optimization, 1998. 13(4): 455-492.
[^2]:  A. Forrester and A. Keane. Engineering design via surrogate modelling: a practical guide. 2008, John Wiley & Sons.
[^3]:  K. Deb. An efficient constraint handling method for genetic algorithms. Computer Methods in Applied Mechanics and Engineering, 2000. 186(2): 311-338.
[^4]:  D. Ginsbourger, R. Le Riche, and L. Carraro. Kriging Is Well-Suited to Parallelize Optimization, in Computational Intelligence in Expensive Optimization Problems, Y. Tenne and C.-K. Goh, Editors. 2010, 131-162.
[^5]:  D. Zhan, J. Qian, and Y. Cheng. Pseudo expected improvement criterion for parallel EGO algorithm. Journal of Global Optimization, 2017. 68(3):  641-662.
[^6]:  C. Chevalier, and D. Ginsbourger. Fast computation of the multi-points expected improvement with applications in batch selection, in Learning and Intelligent Optimization, G. Nicosia and P. Pardalos, Editors. 2013, 59-69.
[^7]: D. Zhan, Y. Meng and H. Xing. A fast multi-point expected improvement for parallel expensive optimization. IEEE Transactions on Evolutionary Computation, 2022, doi: 10.1109/TEVC.2022.3168060.
[^8]:  M. Schonlau. Computer experiments and global optimization. 1997, University of Waterloo.
[^9]: J. Qian, Y. Cheng, J. zhang, J. Liu, and D. Zhan. A parallel constrained efficient global optimization algorithm for expensive constrained optimization problems. Engineering Optimization, 2021. 53(2): 300-320.
[^10]: J. Knowles. ParEGO: A hybrid algorithm with on-line landscape approximation for expensive multiobjective optimization problems. IEEE Transactions on Evolutionary Computation, 2006. 10(1): 50-66.
[^11]: D. Zhan, Y. Cheng, and J. Liu, Expected improvement matrix-based infill criteria for expensive multiobjective optimization. IEEE Transactions on Evolutionary Computation, 2017. 21(6): 956-975.

