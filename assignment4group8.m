% recursive_numerical_methods.m
% - Recursive root-finders: bisection, false-position, Newton-Raphson, secant
% - Knapsack: naive recursive and memoized DP
% - Fibonacci: naive recursive, memoized recursive, iterative DP
% - Time comparisons and plots

clear; close all; clc;

%%  world problem : projectile motion root 
v0 = 20; g = 9.81;
f = @(t) v0.*t - 0.5.*g.*t.^2;
df = @(t) v0 - g.*t;        % derivative for Newton-Raphson
t_exact = 2*v0/g;          % analytical impact time

% Root-search interval and tolerances
a = 0; b = 5;
tol = 1e-8;
maxRecDepth = 1000;        % recursion safety

%% Recursive root finders
% Wrapper functions that measure time and call recursive implementations.

% 1) Recursive Bisection
tic; [root_bis, iter_bis] = recursive_bisection(f, a, b, tol, 0, maxRecDepth); time_bis = toc;

% 2) Recursive False Position (Regula-Falsi)
tic; [root_false, iter_false] = recursive_false_position(f, a, b, tol, 0, maxRecDepth); time_false = toc;

% 3) Recursive Newton-Raphson (recursion on iterations)
x0_nr = 2; tic; [root_newton, iter_newton] = recursive_newton(f, df, x0_nr, tol, 0, maxRecDepth); time_newton = toc;

% 4) Recursive Secant
x0_s = 0.1; x1_s = 4.5; tic; [root_secant, iter_secant] = recursive_secant(f, x0_s, x1_s, tol, 0, maxRecDepth); time_secant = toc;

fprintf('Root-finding results (exact t = %.8f):\n', t_exact);
fprintf('Bisection: root=%.10f, iter=%d, time=%.5e\n', root_bis, iter_bis, time_bis);
fprintf('FalsePos: root=%.10f, iter=%d, time=%.5e\n', root_false, iter_false, time_false);
fprintf('Newton:   root=%.10f, iter=%d, time=%.5e\n', root_newton, iter_newton, time_newton);
fprintf('Secant:   root=%.10f, iter=%d, time=%.5e\n\n', root_secant, iter_secant, time_secant);

%% Plot root-finder computation times
methods = {'Bisection','False Position','Newton-Raphson','Secant'};
times = [time_bis, time_false, time_newton, time_secant];

figure('Name','Root-finder times','NumberTitle','off');
bar(categorical(methods), times);
ylabel('Time (s)');
title('Computation times (recursive implementations)');
grid on;

%% Recursive ODE solvers (single-step recursion)
% Solve dy/dt = -2*y with y(0)=1 on [0,2] for comparison with analytic solution y = exp(-2t)
ode_f = @(t,y) -2*y;
t0 = 0; tf = 2; N = 200; h = (tf - t0)/N;
t_vec = t0:h:tf;
y_exact = exp(-2.*t_vec);

% Euler (recursive)
tic; y_euler = zeros(size(t_vec)); y_euler(1) = 1;
y_euler = recursive_euler_step(ode_f, t_vec, y_euler, 1, h); time_euler = toc;

% Heun (recursive)
tic; y_heun = zeros(size(t_vec)); y_heun(1) = 1;
y_heun = recursive_heun_step(ode_f, t_vec, y_heun, 1, h); time_heun = toc;

% RK4 (recursive)
tic; y_rk4 = zeros(size(t_vec)); y_rk4(1) = 1;
y_rk4 = recursive_rk4_step(ode_f, t_vec, y_rk4, 1, h); time_rk4 = toc;

fprintf('ODE solver timings:\nEuler: %.5e, Heun: %.5e, RK4: %.5e\n\n', time_euler, time_heun, time_rk4);

% Plot ODE results
figure('Name','ODE solvers vs exact','NumberTitle','off');
plot(t_vec, y_exact, 'k-', 'LineWidth', 2); hold on;
plot(t_vec, y_euler, '--r');
plot(t_vec, y_heun, '--b');
plot(t_vec, y_rk4, '--g');
legend('Exact', 'Recursive Euler', 'Recursive Heun', 'Recursive RK4','Location','northeast');
xlabel('t'); ylabel('y(t)'); title('Recursive ODE methods vs analytical');
grid on;

% Bar times
figure('Name','ODE times','NumberTitle','off');
bar(categorical({'Euler','Heun','RK4'}), [time_euler,time_heun,time_rk4]);
ylabel('Time (s)'); title('ODE recursive integrators time (N=200)'); grid on;

%% Knapsack: recursive vs memoized DP 
% Example knapsack instance (small-to-medium size for demonstration)
values = [60, 100, 120, 90, 50];  % values
weights = [10, 20, 30, 40, 10];   % weights
W = 50;                           % capacity
n = length(values);

% 1) Naive recursive (note: exponential time) -- measure time for small n only
tic;
[bestVal_rec, ~] = knapsack_recursive(n, W, weights, values);
time_knap_rec = toc;

% 2) Memoized DP (top-down) 
tic;
memo = -ones(n+1, W+1); % initialize with -1 meaning unknown
[bestVal_mem, memo_used] = knapsack_memo(n, W, weights, values, memo);
time_knap_memo = toc;

% 3) Bottom-up DP (iterative) for verification and faster time
tic;
bestVal_dp = knapsack_dp(weights, values, W);
time_knap_dp = toc;

fprintf('Knapsack results (W=%d, n=%d):\n', W, n);
fprintf('Naive recursion: best=%.2f, time=%.5e\n', bestVal_rec, time_knap_rec);
fprintf('Memoized DP:    best=%.2f, time=%.5e\n', bestVal_mem, time_knap_memo);
fprintf('Bottom-up DP:   best=%.2f, time=%.5e\n\n', bestVal_dp, time_knap_dp);

% Plot knapsack times
figure('Name','Knapsack computation times','NumberTitle','off');
bar(categorical({'Recursive','Memoized','Bottom-up DP'}), [time_knap_rec, time_knap_memo, time_knap_dp]);
ylabel('Time (s)'); title('Knapsack methods timing'); grid on;

%%  Fibonacci: naive vs memo vs iterative DP (timing & scaling) 
% We'll measure for a range of n to show scaling
Ns = [5, 10, 15, 20, 25, 30]; % naive recursion becomes slow quickly (avoid >30)
times_fib_naive = zeros(size(Ns));
times_fib_memo = zeros(size(Ns));
times_fib_iter = zeros(size(Ns));

for idx = 1:length(Ns)
    m = Ns(idx);
    % naive recursion (may be slow for m>30)
    tic; fn = fibonacci_recursive(m); times_fib_naive(idx) = toc;
    % memoized recursive
    tic; fm = fibonacci_memo(m); times_fib_memo(idx) = toc;
    % iterative DP
    tic; fi = fibonacci_iter(m); times_fib_iter(idx) = toc;
end

% Display times
disp(table(Ns', times_fib_naive', times_fib_memo', times_fib_iter', ...
    'VariableNames', {'n','Time_Naive','Time_Memo','Time_Iter'}));

% Plot times (log scale for clarity)
figure('Name','Fibonacci computation times','NumberTitle','off');
semilogy(Ns, times_fib_naive, '-o', 'LineWidth', 1.5); hold on;
semilogy(Ns, times_fib_memo, '-o', 'LineWidth', 1.5);
semilogy(Ns, times_fib_iter, '-o', 'LineWidth', 1.5);
legend('Naive recursion','Memoized recursion','Iterative DP','Location','northwest');
xlabel('n'); ylabel('Time (s) [log scale]'); title('Fibonacci: naive vs memo vs iterative DP'); grid on;

%% Save figures 
 %to save figures to PNG files
 saveas(figure(1),'root_finder_times.png');
 saveas(figure(2),'ode_vs_exact.png');
 saveas(figure(3),'ode_times.png');
 saveas(figure(4),'knapsack_times.png');
saveas(figure(5),'fibonacci_times.png');



%%Nested / Local function implementations 

%% Recursive Bisection
function [c, iter] = recursive_bisection(f, a, b, tol, depth, maxDepth)
    if depth > maxDepth
        error('recursive_bisection: maximum recursion depth reached');
    end
    c = (a + b)/2;
    if abs(f(c)) < tol || (b - a)/2 < tol
        iter = depth;
        return;
    end
    if f(a)*f(c) < 0
        [c, iter] = recursive_bisection(f, a, c, tol, depth+1, maxDepth);
    else
        [c, iter] = recursive_bisection(f, c, b, tol, depth+1, maxDepth);
    end
end

%% Recursive False Position (Regula-Falsi)
function [c, iter] = recursive_false_position(f, a, b, tol, depth, maxDepth)
    if depth > maxDepth
        error('recursive_false_position: maximum recursion depth reached');
    end
    fa = f(a); fb = f(b);
    c = (a*fb - b*fa) / (fb - fa);  % regula-falsi formula
    if abs(f(c)) < tol || abs(b-a) < tol
        iter = depth;
        return;
    end
    if fa * f(c) < 0
        [c, iter] = recursive_false_position(f, a, c, tol, depth+1, maxDepth);
    else
        [c, iter] = recursive_false_position(f, c, b, tol, depth+1, maxDepth);
    end
end

%% Recursive Newton-Raphson (recursion by iteration)
function [x_new, iter] = recursive_newton(f, df, x, tol, depth, maxDepth)
    if depth > maxDepth
        error('recursive_newton: maximum recursion depth reached');
    end
    dfx = df(x);
    if dfx == 0
        error('Derivative zero - Newton method fails at x = %g', x);
    end
    x_new = x - f(x)/dfx;
    if abs(x_new - x) < tol
        iter = depth;
        return;
    else
        [x_new, iter] = recursive_newton(f, df, x_new, tol, depth+1, maxDepth);
    end
end

%% Recursive Secant
function [x2, iter] = recursive_secant(f, x0, x1, tol, depth, maxDepth)
    if depth > maxMaxDepthWrapper(maxDepth)
        error('recursive_secant: maximum recursion depth reached');
    end
    % Avoid division by zero
    denom = f(x1)-f(x0);
    if denom == 0
        x2 = x1; iter = depth; return;
    end
    x2 = x1 - f(x1)*(x1 - x0)/denom;
    if abs(x2 - x1) < tol
        iter = depth;
        return;
    else
        [x2, iter] = recursive_secant(f, x1, x2, tol, depth+1, maxMaxDepthWrapper(maxDepth));
    end
end

function d = maxMaxDepthWrapper(maxDepth)
    % small helper to avoid nested function name conflicts and keep semantics
    d = maxDepth;
end

%% Recursive Euler integrator (fills y vector recursively)
function y = recursive_euler_step(f, t, y, idx, h)
    % idx is the current index (1-based). Recursively compute y(idx+1) until end.
    if idx >= length(t)
        return;
    end
    y(idx+1) = y(idx) + h * f(t(idx), y(idx));
    y = recursive_euler_step(f, t, y, idx+1, h);
end

%% Recursive Heun (improved Euler)
function y = recursive_heun_step(f, t, y, idx, h)
    if idx >= length(t)
        return;
    end
    k1 = f(t(idx), y(idx));
    k2 = f(t(idx) + h, y(idx) + h*k1);
    y(idx+1) = y(idx) + h*(k1 + k2)/2;
    y = recursive_heun_step(f, t, y, idx+1, h);
end

%% Recursive RK4
function y = recursive_rk4_step(f, t, y, idx, h)
    if idx >= length(t)
        return;
    end
    k1 = f(t(idx), y(idx));
    k2 = f(t(idx) + h/2, y(idx) + h*k1/2);
    k3 = f(t(idx) + h/2, y(idx) + h*k2/2);
    k4 = f(t(idx) + h, y(idx) + h*k3);
    y(idx+1) = y(idx) + h*(k1 + 2*k2 + 2*k3 + k4)/6;
    y = recursive_rk4_step(f, t, y, idx+1, h);
end

%% Knapsack implementations
% Naive recursive 0/1 knapsack
function [bestVal, choice] = knapsack_recursive(i, W, weights, values)
    % i = number of items considered (1..n)
    if i == 0 || W == 0
        bestVal = 0;
        choice = [];
        return;
    end
    if weights(i) > W
        [bestVal, choice] = knapsack_recursive(i-1, W, weights, values);
        return;
    else
        % Option 1: don't take item i
        [val1, choice1] = knapsack_recursive(i-1, W, weights, values);
        % Option 2: take item i
        [val2, choice2] = knapsack_recursive(i-1, W - weights(i), weights, values);
        val2 = val2 + values(i);
        if val2 > val1
            bestVal = val2;
            choice = [choice2, i];
        else
            bestVal = val1;
            choice = choice1;
        end
    end
end

% Top-down memoized knapsack
function [bestVal, memo] = knapsack_memo(i, W, weights, values, memo)
    % memo is matrix (i+1) x (W+1), -1 means unknown
    if i == 0 || W == 0
        bestVal = 0;
        return;
    end
    if memo(i+1, W+1) ~= -1
        bestVal = memo(i+1, W+1);
        return;
    end
    if weights(i) > W
        [bestVal, memo] = knapsack_memo(i-1, W, weights, values, memo);
    else
        [val1, memo] = knapsack_memo(i-1, W, weights, values, memo);
        [val2, memo] = knapsack_memo(i-1, W - weights(i), weights, values, memo);
        val2 = val2 + values(i);
        if val2 > val1
            bestVal = val2;
        else
            bestVal = val1;
        end
    end
    memo(i+1, W+1) = bestVal;
end

% Bottom-up DP knapsack (iterative)
function bestVal = knapsack_dp(weights, values, W)
    n = length(weights);
    K = zeros(n+1, W+1);
    for i = 1:n
        for w = 0:W
            if weights(i) <= w
                K(i+1, w+1) = max(K(i, w+1), K(i, w-weights(i)+1) + values(i));
            else
                K(i+1, w+1) = K(i, w+1);
            end
        end
    end
    bestVal = K(n+1, W+1);
end

%% Fibonacci implementations
% Naive recursive Fibonacci
function f = fibonacci_recursive(n)
    if n <= 0
        f = 0; return;
    elseif n == 1
        f = 1; return;
    end
    f = fibonacci_recursive(n-1) + fibonacci_recursive(n-2);
end

% Memoized Fibonacci (top-down)
function f = fibonacci_memo(n)
    memo = -ones(1, n+1);
    f = fib_memo_helper(n, memo);
end

function val = fib_memo_helper(n, memo)
    if n <= 0
        val = 0; return;
    elseif n == 1
        val = 1; return;
    end
    if memo(n+1) ~= -1
        val = memo(n+1); return;
    end
    val = fib_memo_helper(n-1, memo) + fib_memo_helper(n-2, memo);
    memo(n+1) = val; % note: MATLAB passes arrays by value, but it's OK for demonstration small n
end

% Iterative (bottom-up) fibonacci
function f = fibonacci_iter(n)
    if n <= 0
        f = 0; return;
    elseif n == 1
        f = 1; return;
    end
    a = 0; b = 1;
    for k = 2:n
        c = a + b;
        a = b; b = c;
    end
    f = b;
end
