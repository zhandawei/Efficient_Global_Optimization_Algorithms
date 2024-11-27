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
fprintf('Pseuso EI on %d-D %s function, iteration: %d, evaluation: %d, current best solution: %f\n',num_vari,fun_name,iteration,evaluation,fmin);
% the iteration
while evaluation < max_evaluation 
    % build the Kriging model
    kriging_model = Kriging_Train(sample_x,sample_y,lower_bound,upper_bound,1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
    infill_x = zeros(num_q,num_vari);
    point_added = [];
    for ii = 1: num_q
        % find the point with the highest pseudo EI value using GA algorithm
        infill_x(ii,:) = Optimizer_GA(@(x)-Infill_PEI(x,kriging_model,fmin,point_added),num_vari,lower_bound,upper_bound,50,100);
        % update point_added
        point_added = infill_x(1:ii,:);
    end
    % evaluating the query points with the real function 
    best_y = feval(fun_name,infill_x);
    % add the new points to design set
    sample_x = [sample_x;infill_x];
    sample_y = [sample_y;best_y];
    % updating some parameters
    evaluation = size(sample_x,1);
    iteration = iteration + 1;
    fmin = min(sample_y);
    % print the current information to the screen
    fprintf('Pseuso EI on %d-D %s function, iteration: %d, evaluation: %d, current best solution: %f\n',num_vari,fun_name,iteration,evaluation,fmin);
end




