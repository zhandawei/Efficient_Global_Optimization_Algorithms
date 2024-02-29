# A collection of Efficient Global Optimization (EGO) algorithms
## Table of Contents
* [Requirements](#Requirements)
* [Standard EGO](#Standard-Efficient-Global-Optimization-algorithm)
* [Parallel EGO](#Parallel-Efficient-Global-Optimization-algorithms)
* [Constrained EGO](#Constrained-Efficient-Global-Optimization-algorithms)
* [Multiobjecitve EGO](#Multiobjective-Efficient-Global-Optimization-algorithms)
* [To be included](#To-be-included)
* [References](#References)


## Requirements
Windows system.
MATLAB 2016b and above.

## Standard Efficient Global Optimization algorithm
1. **The standard EGO algorithm** (*EGO_EI.m*) [^1]. For the Kriging modeling, the Gaussian correlation function is used as the corrlation function and the constant mean is used as the trend function. I refered some codes in the book *Engineering design via surrogate modelling: a practical guide* [^2] for the Kriging model. The MATLAB **fmincon** function is used for maximizing the likehihood function to get the estimated hyperparameters when training the Kriging model. The expected improvement function is maximized by a real-coded genetic algorithm [^3].

## Parallel Efficient Global Optimization algorithms
1. **The Kriging Believer approach** (*EGO_KB.m*) [^4].
2. **The Constant Liar approach** (*EGO_CL.m*) [^4]. 
3. **The Peseudo Expected Improvement** (*EGO_Pseudo_EI.m*) [^5] .
4. **The Multipoint Expected Improvement** (*EGO_qEI.m*) [^6]. The *qEI* function is coded according the R code in [^7].
5. **The Fast Multipoint Expected Improvement** (*EGO_FqEI.m*) [^8].

## Constrained Efficient Global Optimization algorithms
1. **The Constrained Expected Improvemment** (*EGO_Constrained_EI.m*) [^9].
2. **The Pseudo Constrained Expected Improvement** (*EGO_Pseudo_Constrained_EI.m*) [^10].


## Multiobjective Efficient Global Optimization algorithms
1. **The ParEGO (Pareto EGO)** (*ParEGO.m*) [^11].
2. **The Expected Improvement Matrix** (*EGO_EIM_Euclidean.m*,*EGO_EIM_Hypervolume.m*,*EGO_EIM_Maximin.m*) [^12].


## TO be Included
1. **Multi-EGO**
2. **MOEA/D-EGO**
3. **EGO-EHVI**




## References
[^1]: D. R. Jones, M. Schonlau, and W. J. Welch. Efficient global optimization of expensive black-box functions. Journal of Global Optimization, 1998. 13(4): 455-492.
[^2]:  A. Forrester and A. Keane. Engineering design via surrogate modelling: a practical guide. 2008, John Wiley & Sons.
[^3]:  K. Deb. An efficient constraint handling method for genetic algorithms. Computer Methods in Applied Mechanics and Engineering, 2000. 186(2): 311-338.
[^4]:  D. Ginsbourger, R. Le Riche, and L. Carraro. Kriging Is Well-Suited to Parallelize Optimization, in Computational Intelligence in Expensive Optimization Problems, Y. Tenne and C.-K. Goh, Editors. 2010, 131-162.
[^5]:  D. Zhan, J. Qian, and Y. Cheng. Pseudo expected improvement criterion for parallel EGO algorithm. Journal of Global Optimization, 2017. 68(3):  641-662.
[^6]:  C. Chevalier, and D. Ginsbourger. Fast computation of the multi-points expected improvement with applications in batch selection, in Learning and Intelligent Optimization, G. Nicosia and P. Pardalos, Editors. 2013, 59-69.
[^7]: O. Roustant, D. Ginsbourger, and Y. Deville. DiceKriging, DiceOptim: Two R Packages for the Analysis of Computer Experiments by Kriging-Based Metamodeling and Optimization. Journal of Statistical Software, 2012. 51(1): 1-55.
[^8]: D. Zhan, Y. Meng and H. Xing. A fast multi-point expected improvement for parallel expensive optimization. IEEE Transactions on Evolutionary Computation, 2022, doi: 10.1109/TEVC.2022.3168060.
[^9]:  M. Schonlau. Computer experiments and global optimization. 1997, University of Waterloo.
[^10]: J. Qian, Y. Cheng, J. zhang, J. Liu, and D. Zhan. A parallel constrained efficient global optimization algorithm for expensive constrained optimization problems. Engineering Optimization, 2021. 53(2): 300-320.
[^11]: J. Knowles. ParEGO: A hybrid algorithm with on-line landscape approximation for expensive multiobjective optimization problems. IEEE Transactions on Evolutionary Computation, 2006. 10(1): 50-66.
[^12]: D. Zhan, Y. Cheng, and J. Liu, Expected improvement matrix-based infill criteria for expensive multiobjective optimization. IEEE Transactions on Evolutionary Computation, 2017. 21(6): 956-975.

