classdef (Abstract) NumericalProblem < handle
    % NUMERICALPROBLEM Abstract base class for numerical methods.
    
    properties (Access = protected)
        % Encapsulated properties for the problem definition
        func;   % Function handle to the mathematical function
        params; % Optional parameters for the function
    end
    
    methods
        function obj = NumericalProblem(func, params)
            % Constructor to initialize the function and parameters
            obj.func = func;
            obj.params = params;
        end
    end
    
    methods (Abstract)
        % Abstract method for solving the problem. Subclasses must implement this.
        result = solve(obj);
    end
    
    methods (Access = protected)
        % Encapsulated helper method to validate the function handle
        function validateFunction(obj)
            if ~isa(obj.func, 'function_handle')
                error('The function provided must be a function handle.');
            end
        end
    end
end
