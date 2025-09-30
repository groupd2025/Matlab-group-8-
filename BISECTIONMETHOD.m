%%BISECTION METHOD 
% Ingredients
f=@(x) x^3+3*x-5
a=1;
b=2
n=3;
tic;

SOLUTION

if f(a)*f(b)<0
    for i=1:n
        c=(a+b)/2
        (b-a)<0.01
        b=c
    if f(b)*f(c)<0
        a=c
    end
end
else
    disp('no root between given brackets')
end
I=c
timetaken=toc