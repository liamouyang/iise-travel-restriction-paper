clear
clc

%% Input Configruation
for RRR = 1:8 % Number of cties
    clearvars -except RRR;
    rho = 0.2;
    L = 5; % Number of lockdown levels
    
    Ratio = 0.1; % Basic lockdown budget ratio
    Run_Time_Max = 1800; % main iteration maximum runtime
    Main_Gap = 0.0005; % terminate main iteration

    store = open("file_name_list.mat");
    file_name_list = store.file_name_list;
    open_name = file_name_list{RRR};
    store = open(open_name);

    N = store.N;
    b = store.b;
    c = store.c;
    s = store.s;
    MM = store.MM;
    RHSM = store.RHSM;
    MY = store.MY;
    MX = store.MX;
    RHSS = store.RHSS;
    Grade = store.Grade;
    ProbS = store.ProbS;
    V = store.V;
    Alpha = N; % Hamming distance controller

    %% Initial Solution and Upper Bound
    [x_initial, fei_UB] = Initial_Solution(N, L, V, b, c, MM, RHSM, MY, MX, RHSS, Grade, rho, ProbS);
    % vx_floor: initial feasible solution
    % fei_UB: upper bounds on subgraph variables
    
    %% Define Master Variables
    x = binvar(N, N, L, 'full'); % master variaiable is {0,1}
    x_vector = reshape(x, [N*N*L, 1]); % stretch x to vector
    fei = sdpvar(V, 1); % subgraph variable
    
    %% Initial Hamming Cut
    xx = x_initial; % Assign integer solution by variable rounding
    position0 = find(xx==0);
    position1 = find(xx);
    Hamming = [1 <= sum(x(position0)) + length(position1) - sum(x(position1))];
    
    %% Initial Benders Cut
    Benders = [];
    b_vector = reshape(b, [N*N*L 1]);
    for v = 1:V
        [val_sub, pai, time_sub] = LP_Subproblem(N, L, xx, RHSS{v}, MY, MX, c, s, b);
        Benders = [Benders, pai'*RHSS{v} + (b_vector' - pai'*(MY+MX))*x_vector];
    end
    
    %% Initialziation
    Time_Master = 0;
    Time_Sub_LP = 0;
    Time_Sub_Int = 0;
    RUB = inf;
    RLB = -inf;
    RGap = inf;
    IT = 0; % iteration number
    
    %% Solver Settings
    ops_main = sdpsettings('solver', 'gurobi');
    ops_main.savesolveroutput = 1;
    ops_main.gurobi.TimeLimit = Run_Time_Max;
    ops_main.gurobi.MIPGap = 0.0001;
    
    %% Main Iteration
    % Building Model
    obj_master = rho*sum(b.*x,'all') + (1-rho)*ProbS'*fei; % Benders master problem objective funciotn
    MC = [0<= fei <= fei_UB, MM*x_vector <= RHSM]; % ture master constraint, wont change latter
    
    while Time_Master <= Run_Time_Max && (RGap > Main_Gap || RGap < 0)
        MC_Expand = [MC, Benders, Hamming]; % integerate Benders and Hamming cuts    
        sol_master = optimize(MC_Expand, -obj_master, ops_main);    
        Time_Master = Time_Master + sol_master.solveroutput.result.runtime; % runtime of master problem
        RUB = value(obj_master); % relaxed upper bound
            
        %% Compute Hanmmig Distance Cuts
        xx = value(x); % get master solution
        position0 = find(xx==0);
        position1 = find(xx);
        Hamming = [Hamming, 1 <= sum(x(position0)) + length(position1) - sum(x(position1))];
        
        %% Compute Benders Optimality Cuts
        for v = 1:V
            [LP_Sub_Val(v), pai, time_sub] = LP_Subproblem(N, L, xx, RHSS{v}, MY, MX, c, s, b);
            Benders = [Benders, fei(v) <= pai'*RHSS{v} + (reshape(b, [N*N*L 1])' - pai'*(MY+MX))*x_vector];
            Time_Sub_LP = Time_Sub_LP + time_sub;
        end
        Time_Master = Time_Master + Time_Sub_LP;
    
        %% Compute Relaxed Gap
        RLB = rho*sum(b.*xx,"all") + (1-rho)*ProbS'*LP_Sub_Val'; % Relaxed lower bound
        RGap = (RUB - RLB)/RUB;
    
        IT = IT + 1; % iteration adding 1
    end
    
    %% Phase III Termination
    y_final = cell(V, 1); % store integer subproblem solution
    MC_Final = [MC, Benders]; % final constraint
    
    % Solver Settings
    ops_final = sdpsettings('solver', 'gurobi');
    ops_final.gurobi.TimeLimit = Run_Time_Max;
    ops_final.gurobi.MIPGap = 0.0001;
    ops_final.savesolveroutput = 1;
    
    % Final Master solution
    sol_final = optimize(MC_Final, -obj_master, ops_final);
    x = value(x);
    TrueUB = value(obj_master);
    Time_Master = Time_Master + sol_final.solveroutput.result.runtime;
    
    for v = 1:V
        [Int_Sub_Val(v), time_sub, y_final{v}] = Int_Subproblem(N, L, x, RHSS{v}, MY, MX, c, s, b);
        Time_Sub_Int = Time_Sub_Int + time_sub;
    end
    Time_Master = Time_Master + Time_Sub_Int;
    
    % Compute true gap
    TrueLB =  rho*sum(b.*x, 'all') + (1-rho)*ProbS'*Int_Sub_Val';
    TrueGap = (TrueUB - TrueLB)/TrueUB;
    
    filename = strcat('Results_N', num2str(N), '_L', num2str(L), '_V', num2str(V), '_Ratio', num2str(Ratio), '_Alpha', num2str(Alpha), '_Rho', num2str(rho), '.mat');
    save(filename);

end