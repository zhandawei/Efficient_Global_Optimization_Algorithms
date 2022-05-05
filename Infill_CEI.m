function obj = Infill_CEI(x, kriging_obj,kriging_con,fmin)
% the kriging prediction and varince
[u,s] = Kriging_Predictor(x,kriging_obj);
% the EI value
EI = (fmin-u).*normcdf((fmin-u)./s)+s.*normpdf((fmin-u)./s);
% the number of constraints
num_con = length(kriging_con);
% the kriging prediction and varince
u_g = zeros(size(x,1), num_con);
s_g = zeros(size(x,1), num_con);
for ii = 1: num_con
    [u_g(:, ii), s_g(:, ii)] = Kriging_Predictor(x, kriging_con{ii});
end
% the PoF value
PoF = prod(normpdf((0-u_g)./s_g),2);
CEI = EI.*PoF;
obj = CEI;
end
