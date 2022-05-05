function obj = Infill_EIM_Euclidean(x,kriging_obj,non_dominated_front)
% the input parameters
f = non_dominated_front;
% number of non-dominated points
% number of objectives
[num_pareto,num_obj] = size(f);
% the kriging prediction and varince
u = zeros(1,num_obj);
s = zeros(1,num_obj);
for ii = 1 : num_obj
    [u(:, ii),s(:, ii)] = Kriging_Predictor(x, kriging_obj{ii});
end
u_matrix = repmat(u,num_pareto,1);
s_matrix = repmat(s,num_pareto,1);
EIM = (f - u_matrix).*normcdf((f - u_matrix)./s_matrix) + s_matrix.*normpdf((f - u_matrix)./s_matrix);
y = min(sqrt(sum(EIM.^2,2)));
obj = y;
end
