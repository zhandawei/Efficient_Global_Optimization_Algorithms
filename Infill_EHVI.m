function obj= Infill_EHVI(x,kriging_obj,f)
%-----------------------------------------------------
% number of objectives
num_obj = size(f,2);
r = 1.1*ones(1,num_obj);
% the kriging prediction and varince
u = zeros(1,num_obj);
s = zeros(1,num_obj);
for ii = 1:num_obj
    [u(:, ii),s(:, ii)] = Kriging_Predictor(x,kriging_obj{ii});
end
% EHVI calculated using Monte Carlo simulation
num_simluation_point = 50;
num_simultion_HV = 50;
hvi = zeros(num_simluation_point,1);
rand_sample = mvnrnd(u,diag(s.^2),num_simluation_point);
for ii = 1:num_simluation_point
    if any(all(f <= rand_sample(ii,:),2)) % dominated solution
        hvi(ii) = 0;
    else % non dominated solution
        new_front = [f;rand_sample(ii,:)];
        upper_bound = r;
        lower_bound = min(new_front);
        sim_point = rand(num_simultion_HV,num_obj).*(upper_bound-lower_bound)+lower_bound;
        num_front_point = size(new_front,1);
        simulated_point_matrix = repelem(sim_point,num_front_point,1);
        new_front_matrix = repmat(new_front,num_simultion_HV,1);
        is_dominated = reshape(all(new_front_matrix <= simulated_point_matrix,2),num_front_point,num_simultion_HV)';
        num_improvement = sum(any(is_dominated,2)) - sum(any(is_dominated(:,1:size(f,1)),2));
        hvi(ii) = prod(upper_bound-lower_bound)*num_improvement/num_simultion_HV;
    end
end
obj = mean(hvi);

