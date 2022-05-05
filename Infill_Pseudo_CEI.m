function obj = Infill_Pseudo_CEI(x,kriging_obj,kriging_con,fmin,point_added)
% the kriging prediction and varince
[u,s] = Kriging_Predictor(x,kriging_obj);
% the EI value
EI=(fmin-u).*normcdf((fmin-u)./s)+s.*normpdf((fmin-u)./s);
% the number of constraints
num_con = length(kriging_con);
% the kriging prediction and varince
u_g = zeros(size(x,1), num_con);
s_g = zeros(size(x,1), num_con);
for ii = 1: num_con
    [u_g(:, ii), s_g(:, ii)] = Kriging_Predictor(x, kriging_con{ii});
end
% the PoF value
PoF = prod(normcdf((0-u_g)./s_g), 2);
CEI = EI.*PoF;
lower_bound = kriging_obj.lower_bound;
upper_bound = kriging_obj.upper_bound;
theta = kriging_obj.theta;
% if this is the first infill point
if ~isempty(point_added)
    % the scaling of x is the same for different objectives
    x1 = (x - lower_bound)./(upper_bound- lower_bound);
    x2 = (point_added - lower_bound)./(upper_bound - lower_bound);
    % calculate the correlation
    temp1 = sum(x1.^2.*theta,2)*ones(1,size(x2,1));
    temp2 = sum(x2.^2.*theta,2)*ones(1,size(x1,1));
    correlation = exp(-(temp1 + temp2'-2.*(x1.*theta)*x2'));
    % the Pseudo EI matrix
    CEI = CEI.*prod(1-correlation,2);
end
% the objective is maximized
obj = CEI;
end
