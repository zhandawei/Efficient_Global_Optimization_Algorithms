% The standard EGO algorithm introduced by Jones et al.(1998) [1] for
% solving single-objective bound-constrained expensive optimization
% problems. The EI criterion is maximized by genetic algorithm [2].
% [1] D. R. Jones, M. Schonlau, W. J. Welch. Efficient global optimization
% of expensive black-box functions. Journal of Global Optimization, 1998, 
% 13(4): 455-492.
% [2] K. Deb. An efficient constraint handling method for genetic
% algorithms. Computer Methods in Applied Mechanics and Engineering, 2002,
% 186(2): 311-338.
% This program is free software; you can redistribute it and/or
% modify it. This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
%--------------------------------------------------------------------------
clearvars; close all;
% setting of the problem
fun_name = 'Fun_Rosenbrock';
num_vari = 10;
lower_bound = -2.048*ones(1,num_vari);
upper_bound = 2.048*ones(1,num_vari);
% the number of initial design points
num_initial = 20;
% maximum number of evaluations
max_evaluation = 100;
% initial design points using Latin hypercube sampling method
sample_x = lower_bound + (upper_bound-lower_bound).*lhsdesign(num_initial,num_vari,'criterion','maximin','iterations',1000);
sample_y = feval(fun_name,sample_x);
% the current iteration and evaluation
evaluation = size(sample_x,1);
iteration = 0;
% the current best solution
fmin = min(sample_y);
% print the current information to the screen
fprintf('EGO-EI on %d-D %s function, iteration: %d, evaluation: %d, current best solution: %f\n',num_vari,fun_name,iteration,evaluation,fmin);
% the iteration
while evaluation <  max_evaluation
    % build (or rebuild) the initial Kriging model
    kriging_model = Kriging_Train(sample_x,sample_y,lower_bound,upper_bound,1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
    % find the point with the highest EI value using GA algorithm
    best_x = Optimizer_GA(@(x)-Infill_EI(x,kriging_model,fmin),num_vari,lower_bound,upper_bound,50,100);
    % evaluating the candidate with the real function
    best_y = feval(fun_name,best_x);
    % add the new point to design set
    sample_x = [sample_x;best_x];
    sample_y = [sample_y;best_y];
    % update some parameters
    evaluation = size(sample_x,1);
    iteration = iteration + 1;
    fmin = min(sample_y);
    % print the current information to the screen
    fprintf('EGO-EI on %d-D %s function, iteration: %d, evaluation: %d, current best solution: %f\n',num_vari,fun_name,iteration,evaluation,fmin);
end



