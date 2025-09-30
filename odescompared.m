clc; clear; close all;

f = @(t,y) -2*y;
t0=0; tf=2; y0=1;
N=20; h=(tf-t0)/N;
t=t0:h:tf;

%% Euler
tic;
y_euler=zeros(size(t)); y_euler(1)=y0;
for i=1:N
    y_euler(i+1)=y_euler(i)+h*f(t(i),y_euler(i));
end
time_euler=toc;

%% Heun (Improved Euler)
tic;
y_heun=zeros(size(t)); y_heun(1)=y0;
for i=1:N
    k1=f(t(i),y_heun(i));
    k2=f(t(i)+h,y_heun(i)+h*k1);
    y_heun(i+1)=y_heun(i)+h*(k1+k2)/2;
end
time_heun=toc;

%% RK4
tic;
y_rk4=zeros(size(t)); y_rk4(1)=y0;
for i=1:N
    k1=f(t(i),y_rk4(i));
    k2=f(t(i)+h/2,y_rk4(i)+h*k1/2);
    k3=f(t(i)+h/2,y_rk4(i)+h*k2/2);
    k4=f(t(i)+h,y_rk4(i)+h*k3);
    y_rk4(i+1)=y_rk4(i)+h*(k1+2*k2+2*k3+k4)/6;
end
time_rk4=toc;

%% Analytical Solution
y_exact=exp(-2*t);

%% Plot results
figure;
plot(t,y_exact,'k','LineWidth',2); hold on;
plot(t,y_euler,'--r','LineWidth',1.5);
plot(t,y_heun,'--b','LineWidth',1.5);
plot(t,y_rk4,'--g','LineWidth',1.5);
legend('Exact','Euler','Heun','RK4');
xlabel('t'); ylabel('y(t)');
title('Numerical vs Analytical Solution');
grid on;

%% Display times
fprintf('Computation times:\nEuler: %.6e s\nHeun: %.6e s\nRK4: %.6e s\n', ...
    time_euler,time_heun,time_rk4);
