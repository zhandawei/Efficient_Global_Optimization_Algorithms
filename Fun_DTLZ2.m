function y = Fun_DTLZ2(x,m)
% 0<xi<1   i=1,2...6
% x is a Vector

g=sum((x(:,m:end)-0.5).^2,2);

y(:,1)=(1+g).*prod(cos(0.5*pi*x(:,1:m-1)),2);
if m>2
    for ii=2:m-1
        y(:,ii)=(1+g).*prod(cos(0.5*pi*x(:,1:m-ii)),2).*sin(0.5*pi*x(:,m-ii+1));
    end
end
y(:,m)=(1+g).*sin(0.5*pi*x(:,1));


end