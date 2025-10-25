classdef DifferentialSolver < NumericalProblem
    % DIFFERENTIALSOLVER Solves ordinary differential equations using Euler's method.
    
    properties (Access = private)
        % Private properties specific to differential problems
        tspan; % [t_start, t_end]
        y0;    % Initial condition
        h;     % Step size
    end
    
    methods
        function obj = DifferentialSolver(func, tspan, y0, h)
            % Constructor for a differential problem
            obj = obj@NumericalProblem(func, []); % Call superclass constructor
            obj.tspan = tspan;
            obj.y0 = y0;
            obj.h = h;
            obj.validateFunction(); % Use inherited validation method
        end
        
        function [t_out, y_out] = solve(obj)
            % Implement the Euler's method
            t_start = obj.tspan(1);
            t_end = obj.tspan(2);
            t_out = t_start:obj.h:t_end;
            n_steps = length(t_out);
            y_out = zeros(size(t_out));
            y_out(1) = obj.y0;
            
            for i = 1:(n_steps - 1)
                y_out(i+1) = y_out(i) + obj.h * obj.func(t_out(i), y_out(i));
            end
            
            % Plot the result (optional for visual testing)
            figure;
            plot(t_out, y_out, '-o');
            title('Solution of ODE using Euler''s Method');
            xlabel('t');
            ylabel('y');
        end
    end
end
