% The parallel Efficient Global Optimization (EGO) algorithm [1] using the
% pseudo Expected Improvement criterion (PEI) [2] for solving single-objective 
% unconstrained expensive optimization problems. The PEI criterion is maximized 
% by a geneticalgorithm [3].
% [1] Jones, D.R., Schonlau, M., Welch, W.J.: Efficient global optimization of
%     expensive black-box functions. Journal of Global Optimization 13(4),
%     455-492 (1998).
% [2] D. Zhan, J. Qian, Y. Cheng, Pseudo expected improvement criterion for
% parallel EGO algorithm. Journal of Global Optimization, 2017, 68(3):641-662.
% [3] K. Deb. An efficient constraint handling method for genetic
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
% the number of points selected in each iteration
num_q = 2;
% initial design points using Latin hypercube sampling method
sample_x = repmat(lower_bound,num_initial,1) + repmat(upper_bound-lower_bound,num_initial,1).*lhsdesign(num_initial,num_vari,'criterion','maximin','iterations',1000);
sample_y = feval(fun_name,sample_x);
% the current iteration and evaluation
evaluation = size(sample_x,1);
iteration = 0;
% the current best solution
fmin = min(sample_y);
% print the current information to the screen
fprintf('EGO-PEI on %d-D %s function, iteration: %d, evaluation: %d, current best solution: %f\n',num_vari,fun_name,iteration,evaluation,fmin);
% the iteration
while evaluation < max_evaluation 
    % build (or rebuild) the initial Kriging model
    kriging_model = Kriging_Train(sample_x,sample_y,lower_bound,upper_bound,1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
    % initialize the candidate points and other parameters
    num_k = min(num_q,max_evaluation-evaluation);
    best_x = zeros(num_k,num_vari);
    point_added = [];
    % find the candidates based on pseudo EI criterion
    for ii = 1: num_k
        % find the point with the highest pseudo EI value using GA algorithm
        best_x(ii,:) = Optimizer_GA(@(x)-Infill_Pseudo_EI(x,kriging_model,fmin,point_added),num_vari,lower_bound,upper_bound,50,100);
        % update point_added
        point_added = best_x(1:ii,:);
    end
    % evaluating the candidate with the real function
    best_y = feval(fun_name,best_x);
    % add the new point to design set
    sample_x = [sample_x;best_x];
    sample_y = [sample_y;best_y];
    % updating some parameters
    evaluation = size(sample_x,1);
    iteration = iteration + 1;
    fmin = min(sample_y);
    % print the current information to the screen
    fprintf('EGO-PEI on %d-D %s function, iteration: %d, evaluation: %d, current best solution: %f\n',num_vari,fun_name,iteration,evaluation,fmin);
end




