%% LP Relaxation of Recourse Problem
% CS: RHS vector of subproblems (differ in differ subproblem)
% xx: master solution as input
% val_sub: objective value of subproblem
% pai: dual variables of subproblem
% MY: coefficient matrix of recourse variable y
% MX: coefficient matrix of master decison xx
% c: constraint coefficient matrix in subproblem (same in master)
% s: cost saving mastrix (same in all subproblem)
% b: objective vector or benefit vector (same in all sub and master)

function [val_sub, pai, time_sub] = LP_Subproblem(N, L,  x, CS, MY, MX, c, s, b)

y = sdpvar(N, N, L, 'full'); % y is continuous
obj = sum(b.*y, 'all');
yy = reshape(y, [N*N*L, 1]);
xx = reshape(x, [N*N*L, 1]);
F = [MY*yy <= CS - MX*xx];

% Solver setting
ops = sdpsettings('solver', 'gurobi');
ops.savesolverinput = 1;
ops.savesolveroutput = 1;
ops.gurobi.Method = 1; % use dual simplex to solve LP

warning off; % turn off warning

sol = optimize(F, -obj, ops);

pai = dual(F); % dual variables
val_sub = -sol.solveroutput.result.objval; % dual objective value
time_sub = sol.solveroutput.result.runtime; % runtime

end