%%SIMPSON METHOD
%CREATING AN ANNONYMOUS FUNCTION
f=@(x) x*sin(x)
a=0;
b=pi/2;
n=30;
h=(b-a)/n;
s=f(a)+f(b);
% LOOPING 
for i=1:2:n-2
    s=s+4*f(a+i*h)
end
for i=2:2:n-2
    s=s+2*f(a+i*h)
end
I=s*h/3