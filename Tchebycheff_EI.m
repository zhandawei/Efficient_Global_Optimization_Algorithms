function EI = Tchebycheff_EI(y,s,lamda,z,sample_y)

m = size(y,2);
new_y = lamda(:,1).*(y(1)-z(1));
new_s = lamda(:,1).*s(1);
for ii = 2:m
    mu1 = new_y;
    sigma1_square = new_s.^2;
    mu2 = lamda(:,ii).*(y(ii)-z(ii));
    sigma2_square = lamda(:,ii).^2.*s(ii)^2;
    tao = sqrt(sigma1_square+sigma2_square);
    aphla = (mu1-mu2)./tao;
    new_y = mu1.*gausscdf(aphla) + mu2.*gausscdf(-aphla) + tao.*gausspdf(aphla);
    new_sigma2 = (mu1.^2+sigma1_square).*gausscdf(aphla) + (mu2.^2+sigma2_square).*gausscdf(-aphla) + (mu1+mu2).*gausspdf(aphla) - new_y.^2;
    new_s = sqrt(max(new_sigma2,0));
end

gmin = min(reshape(max(repelem(lamda,size(sample_y,1),1).*repmat(sample_y-z,size(lamda,1),1),[],2),size(sample_y,1),size(lamda,1)))';
EI = (gmin-new_y).*gausscdf((gmin-new_y)./new_s)+new_s.*gausspdf((gmin-new_y)./new_s);
end


function y = gausscdf(x)
y = 0.5*(1+erf(x/sqrt(2)));
end

function res = gausspdf(x)
res = 1/sqrt(2*pi)*exp(-x.^2/2);
end


