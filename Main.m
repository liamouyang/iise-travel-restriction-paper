%% Policy Budget X2
% Risk_Current
% PB_Current
% PC_Later
% PB_Later
% IR_Current
% IR_Later
% ProbS
% rho
clear
clc

N = 21;
L = 5;
rho = 0.6;
delta = 0.05; % policy cost saving ratio

for kkk = 6:6
    load Static_Parameters.mat
    filename = strcat('iteration_', num2str(kkk), '.mat');
    load(filename);
        
    x = binvar(N, N, L, 'full');
    y = binvar(N, N, L, V, 'full');
    
    % Master Problem Constraints
    C = []; % Constraints
    % Constraint (10)
    for i = 1:N
        for j = 1:N
            C = [C, sum(x, 3) == 1];
        end 
    end
    
    % Constraint (11)
    C = [C, sum(PC.*x, [1 3])' <= PB_Current];
    
    % Subproblem Constraints
    for v = 1:V
        C = [C, sum(y(:, :, :, v), 3) == 1]; % Constraint (14)
        for i = 1:N
            for j = 1:N
                for k = 1:L
                   C = [C, sum(y(i, j, [k:L], v)) >= x(i,j,k)]; % Constraint (15)
                end
            end
        end
        C = [C, sum(PC.*y(:, :, :, v), [1 3])' - sum(delta*PC.*x, [1 3])' <= PB_Later{v}]; % Constraint (16)
    end
          
    %% Solver Setting
    ops = sdpsettings('solver', 'gurobi');
    ops.savesolverinput = 1;
    ops.savesolveroutput = 1; % Save original outputs of Gurobi
    ops.gurobi.TimeLimit = 7200;
    ops.gurobi.TuneTimeLimit = 7200;
    % ops.gurobi.Method = 1; % Dual simplex method for LP

    %% OBJ-IR
    Obj_IR = rho*sum(IR_Current.*x, 'all');
    for v = 1:V
        Obj_IR = Obj_IR + (1-rho)* ProbS(v)*sum(IR_Later{v}.*y(:, :, :, v), 'all');
    end
    
    Sol_IR = optimize(C, -Obj_IR, ops);
    x_IR = value(x);
    y_IR = value(y);

    %% Uniform Policy
    C_Uniform = C;
    for i = 1:N
        for j = 1:N-1
            C_Uniform = [ C_Uniform, x(i, j, :) == x(i, j+1, :)];
        end
    end

    for v = 1:V
        for i = 1:N
            for j = 1:N-1
                C_Uniform = [ C_Uniform, y(i, j, :, v) == y(i, j+1, :, v)];
            end
        end
    end
   
    %% OBJ-Uniform-IR
    Sol_IR_Uniform = optimize(C_Uniform, -Obj_IR, ops);
    x_IR_Uniform = value(x);
    y_IR_Uniform = value(y);

    filename = strcat('Results_Risk_Current_Store{', num2str(kkk), '}.mat');
    save(filename);

    clearvars -except kkk N L rho delta
end