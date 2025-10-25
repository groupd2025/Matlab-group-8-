classdef IntegralSolver < NumericalProblem
    % INTEGRALSOLVER Solves definite integrals using the trapezoidal rule.
    
    properties (Access = private)
        % Private properties specific to integral problems
        a; % Lower bound
        b; % Upper bound
        n; % Number of intervals
    end
    
    methods
        function obj = IntegralSolver(func, a, b, n)
            % Constructor for an integral problem
            obj = obj@NumericalProblem(func, []); % Call superclass constructor
            obj.a = a;
            obj.b = b;
            obj.n = n;
            obj.validateFunction(); % Use inherited validation method
        end
        
        function integral_result = solve(obj)
            % Implement the trapezoidal rule
            h = (obj.b - obj.a) / obj.n;
            x = obj.a:h:obj.b;
            y = obj.func(x);
            
            integral_result = h * (0.5 * y(1) + sum(y(2:end-1)) + 0.5 * y(end));
            
            % Display the result
            fprintf('Integral from %f to %f is: %f\n', obj.a, obj.b, integral_result);
        end
    end
end
