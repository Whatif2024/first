classdef Particle < handle
    properties
        position
        velocity
        cost
        best_position
        best_cost
        infeasablity
    end
    
    methods
        function obj = Particle(lower_bound, upper_bound, problem)
            % Constructor
            if nargin > 0
                dim = length(lower_bound);
                obj.position = lower_bound + rand(1, dim) .* (upper_bound - lower_bound);
                obj.velocity = zeros(1, dim);
                
                % Evaluate particle
                [obj.cost, obj.infeasablity] = problem(obj.position);
                
                % Initialize best
                obj.best_position = obj.position;
                obj.best_cost = obj.cost;
            end
        end
        
        function obj = update(obj, w, c, pm, leader, problem)
            % Update particle velocity and position
            dim = length(obj.position);
            
            % Update velocity
            obj.velocity = w * obj.velocity + ...
                          c(1) * rand(1, dim) .* (obj.best_position - obj.position) + ...
                          c(2) * rand(1, dim) .* (leader.position - obj.position);
            
            % Update position
            obj.position = obj.position + obj.velocity;
            
            % Apply mutation
            if rand < pm
                mutation_mask = rand(1, dim) < 0.1; % 10% chance for each dimension
                obj.position(mutation_mask) = obj.position(mutation_mask) + ...
                    0.1 * randn(1, sum(mutation_mask)) .* obj.position(mutation_mask);
            end
            
            % Evaluate new position
            [obj.cost, obj.infeasablity] = problem(obj.position);
            
            % Update personal best
            if obj.infeasablity == 0 && ...
               (obj.best_cost(1) == -1 || dominates(obj.cost, obj.best_cost))
                obj.best_position = obj.position;
                obj.best_cost = obj.cost;
            end
        end
    end
end

function result = dominates(a, b)
    % Check if solution a dominates solution b (for minimization)
    result = all(a <= b) && any(a < b);
end