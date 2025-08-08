classdef Repository < handle
    properties
        swarm
        max_size
        grid_size
        alpha
        beta
        gamma
        grid
        grid_index
    end
    
    methods
        function obj = Repository(initial_swarm, max_size, grid_size, alpha, beta, gamma)
            % Constructor
            obj.max_size = max_size;
            obj.grid_size = grid_size;
            obj.alpha = alpha;
            obj.beta = beta;
            obj.gamma = gamma;
            
            % Initialize with non-dominated feasible particles
            feasible_particles = initial_swarm([initial_swarm.infeasablity] == 0);
            if ~isempty(feasible_particles)
                obj.swarm = obj.findNonDominated(feasible_particles);
            else
                obj.swarm = [];
            end
            
            % Create grid
            obj.createGrid();
        end
        
        function obj = update(obj, new_swarm)
            % Update repository with new particles
            feasible_particles = new_swarm([new_swarm.infeasablity] == 0);
            
            if ~isempty(feasible_particles)
                % Combine with existing repository
                combined = [obj.swarm, feasible_particles];
                
                % Find non-dominated solutions
                obj.swarm = obj.findNonDominated(combined);
                
                % Maintain repository size
                while length(obj.swarm) > obj.max_size
                    obj.deleteOne();
                end
                
                % Update grid
                obj.createGrid();
            end
        end
        
        function leader = SelectLeader(obj)
            % Select leader from repository using grid-based selection
            if isempty(obj.swarm)
                leader = [];
                return;
            end
            
            % Roulette wheel selection based on grid occupancy
            if length(obj.swarm) == 1
                leader = obj.swarm(1);
            else
                % Calculate selection probabilities (prefer less crowded areas)
                grid_counts = hist(obj.grid_index, 1:obj.grid_size^2);
                probs = 1 ./ (grid_counts(obj.grid_index) + 1);
                probs = probs / sum(probs);
                
                % Roulette wheel selection
                cumsum_probs = cumsum(probs);
                r = rand();
                selected_idx = find(cumsum_probs >= r, 1, 'first');
                leader = obj.swarm(selected_idx);
            end
        end
        
        function non_dominated = findNonDominated(obj, particles)
            % Find non-dominated particles
            n = length(particles);
            dominated = false(1, n);
            
            for i = 1:n
                for j = 1:n
                    if i ~= j && dominates(particles(j).cost, particles(i).cost)
                        dominated(i) = true;
                        break;
                    end
                end
            end
            
            non_dominated = particles(~dominated);
        end
        
        function createGrid(obj)
            % Create grid for repository particles
            if isempty(obj.swarm)
                obj.grid = [];
                obj.grid_index = [];
                return;
            end
            
            % Get costs
            costs = vertcat(obj.swarm.cost);
            n_obj = size(costs, 2);
            
            % Find min and max for each objective
            min_vals = min(costs, [], 1);
            max_vals = max(costs, [], 1);
            
            % Expand bounds
            range_vals = max_vals - min_vals;
            min_vals = min_vals - obj.alpha * range_vals;
            max_vals = max_vals + obj.alpha * range_vals;
            
            % Create grid
            obj.grid = cell(1, n_obj);
            for i = 1:n_obj
                obj.grid{i} = linspace(min_vals(i), max_vals(i), obj.grid_size + 1);
            end
            
            % Assign particles to grid cells
            obj.grid_index = zeros(1, length(obj.swarm));
            for i = 1:length(obj.swarm)
                grid_coords = zeros(1, n_obj);
                for j = 1:n_obj
                    grid_coords(j) = find(obj.swarm(i).cost(j) <= obj.grid{j}, 1, 'first') - 1;
                    grid_coords(j) = max(1, min(grid_coords(j), obj.grid_size));
                end
                
                % Convert to linear index
                if n_obj == 2
                    obj.grid_index(i) = (grid_coords(1) - 1) * obj.grid_size + grid_coords(2);
                else
                    obj.grid_index(i) = 1; % Fallback for other dimensions
                end
            end
        end
        
        function deleteOne(obj)
            % Delete one particle from the most crowded grid cell
            if isempty(obj.swarm)
                return;
            end
            
            % Find most crowded grid cell
            grid_counts = hist(obj.grid_index, 1:obj.grid_size^2);
            [~, max_grid] = max(grid_counts);
            
            % Find particles in that grid cell
            candidates = find(obj.grid_index == max_grid);
            
            if length(candidates) == 1
                delete_idx = candidates(1);
            else
                % Random selection among candidates
                delete_idx = candidates(randi(length(candidates)));
            end
            
            % Remove particle
            obj.swarm(delete_idx) = [];
            obj.grid_index(delete_idx) = [];
        end
    end
end

function result = dominates(a, b)
    % Check if solution a dominates solution b (for minimization)
    result = all(a <= b) && any(a < b);
end