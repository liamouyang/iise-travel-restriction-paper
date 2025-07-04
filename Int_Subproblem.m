%% Integer Recourse Problem
% CS: RHS vector of subproblems (differ in differ subproblem)
% xx: master solution as input
% val_sub: objective value of subproblem
% pai: dual variables of subproblem
% MY: coefficient matrix of recourse variable y
% MX: coefficient matrix of master decison xx
% c: constraint coefficient matrix in subproblem (same in master)
% s: cost saving mastrix (same in all subproblem)
% b: objective vector or benefit vector (same in all sub and master)

function [val_sub, time_sub, y] = Int_Subproblem(N, L, x, CS, MY, MX, c, s, b)

y = binvar(N, N, L, 'full'); % y is {0, 1]
obj_IntS = sum(b.*y, 'all');
yy = reshape(y, [N*N*L, 1]);
xx = reshape(x, [N*N*L, 1]);
F = [MY*yy <= CS - MX*xx];

% Solver setting
ops_IntS = sdpsettings('solver', 'gurobi');
% ops_IntS.savesolverinput = 1;
ops_IntS.savesolveroutput = 1;
% ops_IntS.gurobi.TimeLimit = 600;
ops_IntS.gurobi.MIPGap = 0.0001;
ops_IntS.gurobi.NodeMethod = 1; % dual simplex

warning off; % turn off warning
sol_IntS = optimize(F, -obj_IntS, ops_IntS);

val_sub = value(obj_IntS); % dual objective value
time_sub = sol_IntS.solveroutput.result.runtime; % runtime
y = value(y);

end