% Compare recourse value of robust/nominal solution
clear
clc

%% Input Configruation
for RRR = 1:1 % Number of cties
    clearvars -except RRR;
    store = open('file_name.mat');
    open_name = store.file_name{RRR};
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
    V = store.V;
    y_robust = store.y_final;
    rho = store.rho;
    L = store.L;
    x_robust = store.x;

    %% Define Master Variables
    x = binvar(N, N, L, 'full'); % master variaiable is {0,1}
    x_vector = reshape(x, [N*N*L, 1]); % stretch x to vector
    
    %% Gurobi Settings
    ops = sdpsettings('solver', 'gurobi');
    ops.savesolveroutput = 1;
    ops.gurobi.MIPGap = 0.0001;
    
    %% Nominal Current Lockdown Decision
    
    obj_nominal = sum(b.*x,'all'); % original master problem objective funciotn
    MC = [MM*x_vector <= RHSM]; % ture master constraint, wont change latter
    sol_nominal = optimize(MC, -obj_nominal, ops);
    x_nominal = value(x);
    obj_nominal = value(obj_nominal);

    %% Nominal Recourse Decision
    for v = 1:V
        [val_sub_nominal(v), time_sub, y_nominal{v}] = Int_Subproblem(N, L, x_nominal, RHSS{v}, MY, MX, c, s, b);
    end
    val_sub_nominal = val_sub_nominal';

    %% Robust Recourse Decision
    for v = 1:V
        val_sub_robust(v) = sum(b.*y_robust{v}, 'all');
    end

    %% Benefit of Robustness in Recourse
    benefit_recourse = val_sub_robust' - val_sub_nominal;

    %% Benefit of Total
    benefit_total = rho*sum(b.*x_robust, 'all') + (1-rho)*val_sub_robust - rho*sum(b.*x_nominal, 'all') - (1-rho)*val_sub_nominal;

    filename = strcat('Recourse_N', num2str(N), '_L', num2str(L), '_V', num2str(V),  '_Rho', num2str(rho), '.mat');
    save(filename);

end