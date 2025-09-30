%%NEWTON RAPHSON METHOD

% funtion
f=@(x)x^3-4*x-9;
% differentiation
df=@(x)3*x^2-4;
% other given values
e=10^-6;
x0=0;
n=10;
tic;

%SOLUTION

if df(x0~=0)
    for i=1:n
        x1=x0-f(x0)/df(x0);
        x0=x1;
    end
else
    disp('Newton raphson method failed')
end
timetaken=toc;
