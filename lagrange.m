% Known data points
x_values = [0 1 2 3];
y_values = [1 3 2 5];

% Points to interpolate
x_interp = 0:0.1:3;

% Call the function
y_interp = lagrange_interpolation(x_values, y_values, x_interp);

% Plot results
plot(x_values, y_values, 'ro', 'MarkerSize', 8, 'LineWidth', 2); % Known points
hold on
plot(x_interp, y_interp, 'b-', 'LineWidth', 1.5); % Interpolated curve
grid on
xlabel('x');
ylabel('y');
title('Lagrange Interpolation in MATLAB');
legend('Known Points', 'Interpolated Polynomial');
