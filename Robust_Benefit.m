

%% Input Configruation
for RRR = 1:8 % Number of cties
    
    store = open('file_name.mat');
    open_name = store.file_name{RRR};
    store = open(open_name);

    y_robust = store.y_robust;
    x_nominal = store.x_nominal;
    val_sub_nominal = store.val_sub_nominal;
    rho = store.rho;

    val_total_robust = rho*sum(b.*x_robust{RRR}, 'all') + (1-rho)*val_sub_robust;
    val_total_nominal = rho*sum(b.*x_nominal, 'all') + ()




    %% Benefit of Robustness in Recourse
    benefit_recourse = val_sub_robust' - val_sub_nominal;

    %% Benefit of Total
    benefit_total = rho*sum(b.*x_robust, 'all') + (1-rho)*val_sub_robust - rho*sum(b.*x_nominal, 'all') - (1-rho)*val_sub_nominal;

    filename = strcat('Recourse_N', num2str(N), '_L', num2str(L), '_V', num2str(V),  '_Rho', num2str(rho), '.mat');
    save(filename);

end