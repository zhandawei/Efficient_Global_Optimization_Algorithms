%-------------------------------------------------------------------------

clearvars; close all;
% setting of the problem
fun_name = 'Fun_Rosenbrock';
num_vari = 10;
lower_bound = -2.048*ones(1,num_vari);
upper_bound = 2.048*ones(1,num_vari);
% number of initial design points
num_initial = 20;
% maximum number of evaluations
max_evaluation = 100;
% batch size
num_q = 4;
% initial design points using Latin hypercube sampling method
sample_x = lower_bound + (upper_bound-lower_bound).*lhsdesign(num_initial,num_vari,'criterion','maximin','iterations',1000);
sample_y = feval(fun_name,sample_x);
evaluation = size(sample_x,1);
iteration = 0;
% current best solution
fmin = min(sample_y);
% print the current information to the screen
fprintf('Fast Multi-Point Expected Improvement on %d-D %s function, iteration: %d, evaluation: %d, current best solution: %f\n',num_vari,fun_name,iteration,evaluation,fmin);
% the iteration
while evaluation < max_evaluation
    % build the Kriging model
    kriging_model = Kriging_Train(sample_x,sample_y,lower_bound,upper_bound,1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
    % maximize the FqEI function using GA
    [best_x,max_EI]= Optimizer_GA(@(x)-Infill_FqEI(x,kriging_model,fmin),num_vari*num_q,repmat(lower_bound,1,num_q),repmat(upper_bound,1,num_q),50,100);
    infill_x = reshape(best_x,num_vari,[])';
    % evaluate the query points with the real function
    infill_y = feval(fun_name,infill_x);
    % add the new points to design set
    sample_x = [sample_x;infill_x];
    sample_y = [sample_y;infill_y];
    % update some parameters
    evaluation = size(sample_x,1);
    iteration = iteration + 1;
    fmin = min(sample_y);
    % print the current information to the screen
    fprintf('Fast Multi-Point Expected Improvement on %d-D %s function, iteration: %d, evaluation: %d, current best solution: %f\n',num_vari,fun_name,iteration,evaluation,fmin);
end




