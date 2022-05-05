function obj = Infill_Pseudo_PoF(x,kriging_con,point_added)
% the number of constraints
num_con = length(kriging_con);
% the kriging prediction and varince
u = zeros(size(x,1), num_con);
s = zeros(size(x,1), num_con);
for ii = 1: num_con
    [u(:, ii), s(:, ii)] = Kriging_Predictor(x, kriging_con{ii});
end
% the PoF value
PoF = prod(normcdf((0-u)./s), 2);
lower_bound = kriging_con{1}.lower_bound;
upper_bound = kriging_con{1}.upper_bound;
theta = kriging_con{1}.theta;
% if this is the first infill point
if ~isempty(point_added)
    % the scaling of x
    x1 = (x - lower_bound)./(upper_bound- lower_bound);
    x2 = (point_added - lower_bound)./(upper_bound - lower_bound);
    % calculate the correlation
    temp1 = sum(x1.^2.*theta,2)*ones(1,size(x2,1));
    temp2 = sum(x2.^2.*theta,2)*ones(1,size(x1,1));
    correlation = exp(-(temp1 + temp2'-2.*(x1.*theta)*x2'));
    % the Pseudo EI matrix
    PoF=PoF.*prod(1-correlation,2);
end
% the objective is maximized
obj = PoF;
end
