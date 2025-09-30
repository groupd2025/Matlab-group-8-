%%FIXED POINT ITERATION
%FIXED POINT ITERATION

% annonymous function()
f=@(x) (2^x+2)/5;
x0=0;
e=0.0001;
n=10;
tic;

%5SOLUTION

for i=1:n
    x1=f(x0);
    fprintf('x%d=%4f\n',i,x1)
    if abs(x1-x0)<e
    end
    x0=x1;
end
timetaken=toc;

