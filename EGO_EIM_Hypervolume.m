% The multiobjective EGO algorithm using EIM(expected improvement
%    matrix)-based criteria, which is significant cheaper-to-evaluate than the
%    state-of-the-art multiobjective EI criteria. For detailed description
%    about the EIM criteria, please refer to [1]. The non-dominated sorting method
% by Yi Cao [2] is used to identify the non-dominated fronts from all the design points
% The hypervolume indicators are calculated using the faster algorithm of
%     Nicola Beume et al. (2009) [3]. The EIM criteria are maximized by GA [4] algorithm.
% [1]  D. Zhan, Y. Cheng, J. Liu, Expected Improvement Matrix-based Infill
%      Criteria for Expensive Multiobjective Optimization. IEEE Transactions
%      on Evolutionary Computation, 2017, 21 (6): 956-975.
% [2] http://www.mathworks.com/matlabcentral/fileexchange/17251-
%      pareto-front.
% [3] N. Beume, C.M. Fonseca, M. Lopez-Ibanez, L. Paquete, J. Vahrenhold,
%     On the Complexity of Computing the Hypervolume Indicator, IEEE
%     Transactions on Evolutionary Computation 13(5) (2009) 1075-1082.
% [4] K. Deb. An efficient constraint handling method for genetic
% algorithms. Computer Methods in Applied Mechanics and Engineering, 2002,
% 186(2): 311-338.
% This program is free software; you can redistribute it and/or
% modify it. This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
% -----------------------------------------------------------------------------------------
clearvars;close all;
% settings of the problem
fun_name = 'Fun_DTLZ2';
% number of objectives
num_obj = 2;
% number of design variables
num_vari = 6;
lower_bound = zeros(1,num_vari);
upper_bound = ones(1,num_vari);
% reference point for calculating hypervolume
ref_point = 2.5*ones(1,num_obj);
% number of initial design points
num_initial = 20;
% maximum number of evaluations
max_evaluation = 100;
% the intial design points, points sampled all at once
sample_x = lower_bound + (upper_bound-lower_bound).*lhsdesign(num_initial,num_vari,'criterion','maximin','iterations',1000);
sample_y = feval(fun_name, sample_x, num_obj);
% scale the objectives to [0,1]
sample_y_scaled =(sample_y - min(sample_y))./(max(sample_y)-min(sample_y));
% initialize some parameters
iteration = 0;
evaluation = size(sample_x,1);
kriging_obj = cell(1,num_obj);
% calculate the initial hypervolume values
index = Paretoset(sample_y);
non_dominated_front = sample_y(index,:);
non_dominated_front_scaled = sample_y_scaled(index,:);
hypervolume = Hypervolume(non_dominated_front,ref_point);
% print the hypervolume information
fprintf('EGO-EIM-Hypervolume on %d-D %s function, iteration: %d, evaluation: %d, hypervolume: %f \n',num_vari,fun_name,iteration,evaluation,hypervolume);
% beginning of the iteration
while evaluation < max_evaluation
    % build the initial kriging model for each objective
    for ii = 1:num_obj
        kriging_obj{ii} = Kriging_Train(sample_x,sample_y_scaled(:,ii),lower_bound,upper_bound,1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
    end
    % select updating points using the EIM criteria
    best_x = Optimizer_GA(@(x)-Infill_EIM_Hypervolume(x,kriging_obj,non_dominated_front_scaled),num_vari,lower_bound,upper_bound,50,100);
    best_y = feval(fun_name,best_x, num_obj);
    % add the new points to the design set
    sample_x = [sample_x;best_x];
    sample_y = [sample_y;best_y];
    sample_y_scaled = (sample_y - min(sample_y))./(max(sample_y)-min(sample_y));
    evaluation = evaluation + size(best_x,1);
    iteration = iteration + 1;
    % calculate the hypervolume values
    index = Paretoset(sample_y);
    non_dominated_front = sample_y(index,:);
    non_dominated_front_scaled = sample_y_scaled(index,:);
    hypervolume = Hypervolume(non_dominated_front,ref_point);
    % plot current non-dominated front points
    % print the hypervolume information
    fprintf('EGO-EIM-Hypervolume on %d-D %s function, iteration: %d, evaluation: %d, hypervolume: %f\n',num_vari,fun_name,iteration,evaluation,hypervolume);
end

