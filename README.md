# Matlab implementations of Bayesian Optimization algorithms
Bayesian Optimization (BO) algorithms, also known as Efficient Global Optimization (EGO) algorithms are widely used to solve expensive optimization problems. I try to make a collection of different Bayesian optimization algorithms that we have proposed and implemented during my research. In these implementations, I try keep the codes as simple as possile.


## Table of Contents
* [Requirements](#Requirements)
* [Standard Algorithm](#Standard-algorithm)
* [High-Dimensional Algorithms](#High-Dimensional-Algorithms)
* [Parallel (Batch) Algorithms](#Parallel-Efficient-Global-Optimization-algorithms)
* [Multiobjecitve Algorithms](#Multiobjective-Efficient-Global-Optimization-algorithms)
* [Constrained Algorithms](#Constrained-Efficient-Global-Optimization-algorithms)
* [References](#References)


## Requirements
1. **Windows system.** I have not tested other operating systems, but the codes should also work as Matlab is cross-platformed.
2. **MATLAB 2016b and above**. I used a lot of ```.*``` to multiply vector and matrix. This multiplication uses implicit expansion, which was introduced in MATLAB 2016b.

## Standard Bayesian Optimization algorithm
**The standard BO algorithm**[^1] ```Standard_BO.m``` .
 
For the Kriging modeling, the Gaussian correlation function is used as the corrlation function and the constant mean is used as the trend function. 
I refered some codes in the book *Engineering design via surrogate modelling: a practical guide* [^2] for the Kriging model. 
The MATLAB ```fmincon``` function is used for maximizing the likehihood function to get the estimated hyperparameters when training the Kriging model. 
The expected improvement function is maximized by a real-coded genetic algorithm [^3].


## High-Dimensional Bayesian Optimization Algorithms
1. **The Dropout Approach**[^4] ```HD_Dropout.m``` 



## Parallel Bayesian Optimization Algorithms
1. **The Kriging Believer Approach** ```Parallel_KB.m``` [^4].
   
3. **The Constant Liar Approach** (*Parallel_CL.m*) [^4].
4. 
5. **The Peseudo Expected Improvement** (*Parallel_PEI.m*) [^5] .
6. 
7. **The Multipoint Expected Improvement** (*Parallel_qEI.m*) [^6].
8. The *qEI* function is coded according the R code in [^7].
9. 
10. **The Fast Multipoint Expected Improvement** (*Parallel_FqEI.m*) [^8].
11. 

## Multiobjective Bayesian Optimization Algorithms
1. **The ParEGO (Pareto EGO)** (*Multiobjective_ParEGO.m*) [^9].
2. **The Expected Improvement Matrix** (*Multiobjective_EIM_Euclidean.m*,*Multiobjective_EIM_Hypervolume.m*,*Multiobjective_EIM_Maximin.m*) [^10].
3. **The Expected Hypervolume Improvement** (*Multiobjective_EHVI.m*)[^11]. The *EHVI* criterion is calculated using Monte Carlo approximation.
4. **The MOEA/D-EGO** (*Multiobjective_MOEAD_EGO.m*)[^12]. We use all the samples to train the Kriging models instead of using the fuzzy clusting based modeling method used in the original work[^12].


## Constrained Bayesian Optimization Algorithms
1. **The Constrained Expected Improvemment** (*Constrained_CEI.m*) [^13].
2. **The Pseudo Constrained Expected Improvement** (*Constrained_PCEI.m*) [^14].

## References
[^1]: D. R. Jones, M. Schonlau, and W. J. Welch. Efficient global optimization of expensive black-box functions. Journal of Global Optimization, 1998. 13(4): 455-492.
[^2]:  A. Forrester and A. Keane. Engineering design via surrogate modelling: a practical guide. 2008, John Wiley & Sons.
[^3]:  K. Deb. An efficient constraint handling method for genetic algorithms. Computer Methods in Applied Mechanics and Engineering, 2000. 186(2): 311-338.
[^4]: C. Li, S. Gupta, S. Rana, T. V. Nguyen, S. Venkatesh, and A. Shilton. High dimensional bayesian optimization using dropout, in International Joint Conference on Artificial Intelligence, 2017, 2096-2102.
[^5]:  D. Ginsbourger, R. Le Riche, and L. Carraro. Kriging Is Well-Suited to Parallelize Optimization, in Computational Intelligence in Expensive Optimization Problems, Y. Tenne and C.-K. Goh, Editors. 2010, 131-162.
[^5]:  D. Zhan, J. Qian, and Y. Cheng. Pseudo expected improvement criterion for parallel EGO algorithm. Journal of Global Optimization, 2017. 68(3):  641-662.
[^6]:  C. Chevalier, and D. Ginsbourger. Fast computation of the multi-points expected improvement with applications in batch selection, in Learning and Intelligent Optimization, G. Nicosia and P. Pardalos, Editors. 2013, 59-69.
[^7]: O. Roustant, D. Ginsbourger, and Y. Deville. DiceKriging, DiceOptim: Two R Packages for the Analysis of Computer Experiments by Kriging-Based Metamodeling and Optimization. Journal of Statistical Software, 2012. 51(1): 1-55.
[^8]: D. Zhan, Y. Meng and H. Xing. A fast multi-point expected improvement for parallel expensive optimization. IEEE Transactions on Evolutionary Computation, 2022, doi: 10.1109/TEVC.2022.3168060.
[^9]: J. Knowles. ParEGO: A hybrid algorithm with on-line landscape approximation for expensive multiobjective optimization problems. IEEE Transactions on Evolutionary Computation, 2006. 10(1): 50-66.
[^10]: D. Zhan, Y. Cheng, and J. Liu. Expected improvement matrix-based infill criteria for expensive multiobjective optimization. IEEE Transactions on Evolutionary Computation, 2017. 21(6): 956-975.
[^11]: M. T. M. Emmerich, K. C. Giannakoglou, and B. Naujoks. Single- and multiobjective evolutionary optimization assisted by Gaussian random field metamodels. IEEE Transactions on Evolutionary Computation, 2006, 10(4): 421-439.
[^12]: Q. Zhang, W. Liu, E. Tsang, and B. Virginas. Expensive Multiobjective Optimization by MOEA/D With Gaussian Process Model. IEEE Transactions on Evolutionary Computation, 2010, 14(3): 456-474.
[^13]:  M. Schonlau. Computer experiments and global optimization. 1997, University of Waterloo.
[^14]: J. Qian, Y. Cheng, J. zhang, J. Liu, and D. Zhan. A parallel constrained efficient global optimization algorithm for expensive constrained optimization problems. Engineering Optimization, 2021. 53(2): 300-320.
