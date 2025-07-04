% Phase I Heuristic Initial Cuts

function [vx_floor, fei_UB] = Initial_Solution(N, L, V, b, c, MM, CM, MY, MX, Budget_Sub, Grade, rho, ProbS)
    %% Phase 1.1. Solving LP Relaxation of Orignial Problem
    x = sdpvar(N, N, L, 'full');
    y = sdpvar(N, N, L, V, 'full');
    xxx = reshape(x, [N*N*L 1]);
    yyy = cell(V);
    
    for v = 1:V
        yyy{v} = reshape(y(:, :, :, v), [N*N*L 1]);
    end
    
    obj = rho*sum(b.*x, 'all');
    for v = 1:V
        obj = obj + (1-rho)*ProbS(v)*sum(b.*y(:, :, :, v), 'all');
    end
    
    C_Initial = [0 <= x <=1]; % Add Variable Bounds on x!!!
    C_Initial =[C_Initial, MM*xxx <= CM];

    for v = 1:V
        C_Initial = [C_Initial, MY*yyy{v} + MX*xxx <= Budget_Sub{v}];
    end
    
    %% Solver Setting
    ops = sdpsettings('solver', 'gurobi');
    ops.savesolverinput = 1;
    ops.savesolveroutput = 1;
    ops.gurobi.Method = 1; % use dual simplex to solve LP
    sol = optimize(C_Initial, -obj, ops);
    
    %% Phase 1.2. Initial Subgraph Variable Upper Bound
    for v = 1:V
        fei_UB(v) = sum(b.*value(y(:, :, :, v)), 'all');
    end
    fei_UB = fei_UB';
    
    %% Phase 1.3. Fractional Variable Rounding
    % find index for fractional lockdown level
    vx = value(x);
    level_fra = zeros(N, N);
    for j = 1:N
        for i = 1:N
            k = L;
            is_fractional = 0;
            while is_fractional == 0 && k >= 1
                if vx(i, j, k) > 0 && vx(i, j, k) < 1
                    level_fra(i, j) = k; % recrod the positon of fractional variable
                    is_fractional = 1;
                else
                    k = k - 1;
                end
            end        
        end
    end
    
    % find inbound city list with fractional lockdown level
    city_fra = cell(N, 1);
    for j = 1:N
        city_fra{j} = find(level_fra(:, j));
    end
    
    % Extract Grade of fraction lockdown variable
    for j = 1:N
        for i = 1:N
            if level_fra(i, j) > 0
                true_grade(i, j) = Grade(i, j, level_fra(i, j));
            else
                true_grade(i, j) = 0; 
            end
        end
    end

    % Rounding from city list entry of largest true_grade
    vx_floor = floor(vx); % round down all fractional variables to be 0
    cost_current = sum(c.*vx_floor, [1 3])'; % used cost by vx_floor
    cost_remaining = full(CM(end-N+1:end)) - cost_current; % remainging cost for rounding

    % Ranking grade from largest to smallest and return city index
     for j = 1:N
        [loca, true_position] = ismember(sort(true_grade(:, j), 'descend'), true_grade(:, j));
        % extract true city index of k-largest grade

        k = 1; % from the city of k-largest grade
        tloc = true_position(k);

        while true_grade(tloc, j) > 0
            % true_position{j}(k) is the k-largest grade
            % return the true grade
            % if the grade is 0, break
           
            if c(tloc, j, level_fra(tloc, j)) <= cost_remaining(j) 
            % rounding the fractional variable up does not exceed remaining budget
          
            vx_floor(tloc, j, level_fra(tloc, j)) = 1; 
            % position(k) is the true city index (not grade index)           
            
            cost_remaining(j) = cost_remaining(j) - c(tloc, j, level_fra(tloc, j)); 
            % update the remaining cost            
            end   
            
            % check the city of (k+1)-largest grade
            k = k + 1;
            tloc = true_position(k);
        end
     end

end