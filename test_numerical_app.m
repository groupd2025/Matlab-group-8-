% test_numerical_app.m

% --- Test the DifferentialSolver subclass ---
fprintf('--- Testing DifferentialSolver (Euler''s Method) ---\n');

% Define the ODE: dy/dt = -2t*y
ode_func = @(t,y) -2*t.*y;
tspan_ode = [0, 2];
y0_ode = 1;
h_ode = 0.1;

% Create an object of the DifferentialSolver class
diff_problem = DifferentialSolver(ode_func, tspan_ode, y0_ode, h_ode);

% Solve the differential problem using the polymorphic 'solve' method
[t_sol, y_sol] = diff_problem.solve();

% The analytical solution is y(t) = exp(-t^2)
y_exact = exp(-t_sol.^2);

% Plot the numerical and analytical solutions for comparison
hold on;
plot(t_sol, y_exact, 'r-', 'LineWidth', 2);
legend('Euler Method', 'Exact Solution');
hold off;
fprintf('Differential problem solved. Plot generated.\n\n');

% --- Test the IntegralSolver subclass ---
fprintf('--- Testing IntegralSolver (Trapezoidal Rule) ---\n');

% Define the function to integrate: f(x) = x^2 from 0 to 1
integral_func = @(x) x.^2;
a_int = 0;
b_int = 1;
n_intervals = 100;

% Create an object of the IntegralSolver class
int_problem = IntegralSolver(integral_func, a_int, b_int, n_intervals);

% Solve the integral problem using the polymorphic 'solve' method
integral_result = int_problem.solve();

% The analytical solution is 1/3
fprintf('Expected integral result is 0.333333...\n');
