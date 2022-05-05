function obj = Infill_EI(x,Kriging_model,fmin)
% get the Kriging prediction and variance
[u,s] = Kriging_Predictor(x,Kriging_model);
% calcuate the EI value
EI = (fmin-u).*normcdf((fmin-u)./s)+s.*normpdf((fmin-u)./s);
% this EI needs to be maximized
obj = EI;

end





