% clear
% clc

%% Input Configruation
for N = [30] % Number of cties
    rho = 0.5;
    L = 5; % Number of lockdown levels
    Alpha = N; % Hamming distance controller
    Ratio = 0.1; % Basic lockdown budget ratio
    
    %% Define Master Variables
    x = binvar(N, N, L, 'full'); % master variaiable is {0,1}
    y = binvar(N, N, L, 'full');
    x_vector = reshape(x, [N*N*L, 1]); % stretch x to vector
    y_vector = reshape(y, [N*N*L, 1]);      
        
    %% Main Iteration
    % Building Model
    obj = rho*sum(b.*x,'all') + (1-rho)*sum(b.*y,'all'); 
    MC = [MM*x_vector <= RHSM]; % ture master constraint, wont change latter
    MC = [MC, MY*y_vector + MX*x_vector <= RHSS{1}];

    %% Solver Settings
    ops_main = sdpsettings('solver', 'gurobi');
    ops_main.savesolveroutput = 1;
    ops_main.gurobi.MIPGap = 0.0001;

    sol = optimize(MC, -obj, ops_main);
    x_old_robust = value(x);
    y_old_robust = value(y);
    val_old_robust = value(obj);
    
    filename = strcat('New_Results_N', num2str(N), '_L', num2str(L), '_V', num2str(V), '_Ratio', num2str(Ratio), '_Alpha ', num2str(Alpha), '.mat');
    save(filename);

end