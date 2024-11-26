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
sample_x = lower_bound + (upper_bound-lower_bound).*lhsdesign(num_initial,num_vari,'criterion','maximin','iterations',1000);
sample_y = feval(fun_name,sample_x);
% the current iteration and evaluation
evaluation = size(sample_x,1);
iteration = 0;
% the current best solution
fmin = min(sample_y);
% print the current information to the screen
fprintf('KB on %d-D %s function, iteration: %d, evaluation: %d, current best solution: %f\n',num_vari,fun_name,iteration,evaluation,fmin);
% the iteration
while evaluation < max_evaluation
    % build (or rebuild) the initial Kriging model
    kriging_model = Kriging_Train(sample_x,sample_y,lower_bound,upper_bound,1*ones(1,num_vari),0.001*ones(1,num_vari),1000*ones(1,num_vari));
    % maximize the qEI function
    infill_x = zeros(num_q,num_vari);
    for ii = 1: num_q
        % the Expected Improvement criterion
        [infill_x(ii,:),max_EI] = Optimizer_GA(@(x)-Infill_EI(x,kriging_model,fmin),num_vari,lower_bound,upper_bound,num_vari,100);
        % use kriging prediction as the real objective
        kriging_model = Kriging_Train([sample_x;infill_x(1:ii,:)],[sample_y;Kriging_Predictor(infill_x(1:ii,:),kriging_model)],lower_bound,upper_bound,kriging_model.theta);
    end
    % evaluating the candidates with the real function
    infill_y = feval(fun_name,infill_x);
    % add the new point to design set
    sample_x = [sample_x;infill_x];
    sample_y = [sample_y;infill_y];
    % updating some parameters
    evaluation = size(sample_x,1);
    iteration = iteration + 1;
    fmin = min(sample_y);
    % print the current information to the screen
    fprintf('KB on %d-D %s function, iteration: %d, evaluation: %d, current best solution: %f\n',num_vari,fun_name,iteration,evaluation,fmin);
end




