
rho = 0.8; % REMEMBER Changing RHO!+
L = 5;
x_nominal_all = cell(8, 1);

for RRR = 1:8
    
    %% Obtain traditional robust solution
    N = N_all{RRR};
    MM = MM_all{RRR};
    RHSM = RHSM_all{RRR};
    b = b_all{RRR};
    MX = MX_all{RRR};
    MY = MY_all{RRR};
    RHSS = RHSS_all{RRR};
    RHSS = RHSS{1};

    x = binvar(N, N, L, 'full');
    y = binvar(N, N, L, 'full');

    obj = rho*sum(b.*x, 'all') + (1-rho)*sum(b.*y, 'all');
    F = [MM*reshape(x, [N*N*L 1]) <= RHSM, MX*reshape(x, [N*N*L 1]) + MY*reshape(y, [N*N*L 1]) <= RHSS];
    ops = sdpsettings('solver', 'gurobi');
    sol = optimize(F, -obj, ops);
    x_nominal_all{RRR} = value(x);

    %% Compare with DRO solution
    RHSS = RHSS_all{RRR};
    V = V_all{RRR};
    x_robust = x_robust_all{RRR};
    x_nominal = x_nominal_all{RRR};
    val_robust = zeros(V, 1);
    val_nominal = zeros(V, 1);

    for v = 1:V
        [val_robust(v), time_sub, y_robust] = Int_Subproblem(N, L, x_robust, RHSS{v}, MY, MX, b);
        [val_nominal(v), time_sub, y_nominal] = Int_Subproblem(N, L, x_nominal, RHSS{v}, MY, MX, b);
    end

    val_robust_all{RRR} = val_robust;
    val_nominal_all{RRR} = val_nominal;

    %% Find good instances
%     diff = val_robust - val_nominal;
%     ind_diff = [];
%     for v = 1:V
%         if diff(v) > 0
%             ind_diff = [ind_diff, v];
%         end
%     end

end

filename = 'Resultt_Robust_New_RHO08.mat';
save(filename);