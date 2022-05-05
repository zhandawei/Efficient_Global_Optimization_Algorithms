function obj = Infill_PoF(x, kriging_con)
% the number of constraints
num_con = length(kriging_con);
% the kriging prediction and varince
u = zeros(size(x,1), num_con);
s = zeros(size(x,1), num_con);
for ii = 1: num_con
    [u(:, ii), s(:, ii)] = Kriging_Predictor(x, kriging_con{ii});
end
% the PoF value
PoF = prod(normcdf((0-u)./s),2);
obj = PoF;
end
