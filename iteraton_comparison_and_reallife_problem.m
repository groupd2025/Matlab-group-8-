clc; clear; close all;

% Function (height of projectile)
v0 = 20; g = 9.81;
f = @(t) v0*t - 0.5*g*t.^2;
fprime = @(t) v0 - g*t;   % derivative for Newton-Raphson

% Exact solution (impact time)
t_exact = 2*v0/g;

% Parameters
tol = 1e-6; maxIter = 100;
a=1; b=5;      % interval containing root
x0=2;          % initial guess

methods = {'Bisection','False Position','Newton-Raphson','Secant'};
roots_found = zeros(1,length(methods));
times = zeros(1,length(methods));

%% Bisection Method
tic;
a1=a; b1=b;
for i=1:maxIter
    c=(a1+b1)/2;
    if abs(f(c))<tol, break; end
    if f(a1)*f(c)<0, b1=c; else a1=c; end
end
roots_found(1)=c; times(1)=toc;

%% False Position Method
tic;
a1=a; b1=b;
for i=1:maxIter
    c=(a1*f(b1)-b1*f(a1))/(f(b1)-f(a1));
    if abs(f(c))<tol, break; end
    if f(a1)*f(c)<0, b1=c; else a1=c; end
end
roots_found(2)=c; times(2)=toc;

%% Newton-Raphson Method
tic;
x=x0;
for i=1:maxIter
    x_new = x - f(x)/fprime(x);
    if abs(x_new-x)<tol, break; end
    x=x_new;
end
roots_found(3)=x_new; times(3)=toc;

%% Secant Method
tic;
x0=1; x1=5;
for i=1:maxIter
    x2 = x1 - f(x1)*(x1-x0)/(f(x1)-f(x0));
    if abs(x2-x1)<tol, break; end
    x0=x1; x1=x2;
end
roots_found(4)=x2; times(4)=toc;

%% Display results
fprintf('\nProjectile Motion Root-finding Comparison:\n');
for i=1:length(methods)
    fprintf('%s: root = %.6f, time = %.6e s, error = %.2e\n', ...
        methods{i}, roots_found(i), times(i), abs(roots_found(i)-t_exact));
end
fprintf('Exact solution: %.6f s\n', t_exact);

%% Plot function and roots
t_vals=linspace(0,5,200);
figure;
plot(t_vals,f(t_vals),'b','LineWidth',2); hold on;
yline(0,'k--');
plot(roots_found,zeros(size(roots_found)),'ro','MarkerSize',8,'MarkerFaceColor','r');
plot(t_exact,0,'gs','MarkerSize',10,'MarkerFaceColor','g');
legend('f(t)','y=0','Approx roots','Exact root');
xlabel('t (s)'); ylabel('f(t)'); title('Projectile Motion Root-Finding');
grid on;

%% Bar chart of computation times
figure;
bar(categorical(methods),times);
ylabel('Time (s)'); title('Computation Times of Methods');
