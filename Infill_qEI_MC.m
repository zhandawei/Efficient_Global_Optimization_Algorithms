function EI = Infill_qEI_MC(new_x,model,fmin)
n = 1E5;
% reshape the input
new_x = reshape(new_x,size(model.S,2),[])';
% predictions and errors
[u,~,~,Cov] = predictor(new_x,model);
mu = u';
sigma = Cov;
sample = mvnrnd(mu,sigma,n);
EI  = mean(max(fmin - min(sample,[],2),0));