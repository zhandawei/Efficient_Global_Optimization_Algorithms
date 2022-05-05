function [obj, con] = Fun_Welded_Beam(x)
% ----------------------------------------------------------------------------
% the welded beam design problem
% Rao, S.S. Engineering Optimization. Wiley, New York. 1996.
% fmin = 1.7250
% ----------------------------------------------------------------------------
x1 = x(:, 1); x2 = x(:, 2); x3 = x(:, 3); x4 = x(:, 4);

P = 6000;
L = 14;
E = 30*10^6;
G = 12 * 10^6;
tao_max = 13600;
sigma_max = 30000;
deta_max = 0.25;

M = P*(L + 0.5*x2);
R = sqrt(x2.^2/4 + ((x1+x3)/2).^2);
J = 2*sqrt(2)*x1.*x2.*(x2.^2/12 + ((x1+x3)/2).^2);
tao_1 = P./(sqrt(2)*x1.*x2);
tao_2 = M.*R./J;

tao = sqrt(tao_1.^2 + 2*tao_1.*tao_2.*x2./(2*R) + tao_2.^2);
sigma = 6*P*L./(x4.*x3.^2);
deta = 4*P*L.^3./(E*x3.^3.*x4);
Pc = 4.013*E*sqrt(x3.^2.*x4.^6/36).*(1-x3.*sqrt(E/(4*G))/(2*L))/L^2;

obj = 1.10471*x1.^2.*x2 + 0.04811*x3.*x4.*(x2 + 14);
con(:, 1) = tao./tao_max - 1;
con(:, 2) = sigma./sigma_max - 1;
con(:, 3) = (x1 - x4)/10;
con(:, 4) = 0.10471*x1.^2 +0.04811*x3.*x4.*(14.0+x2) -5;
con(:, 5) = deta./deta_max - 1;
con(:, 6) = 1 - Pc./P;



end