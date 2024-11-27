function f = Rosenbrock(x)
% the Rosenbrock function
% xi = [-2.048,2.048]
f = sum(100*(x(:,2:end) - x(:,1:end-1).^2).^2 + (x(:,1:end-1)-1).^2,2);
end