% clear
% clc
% 
% %% Input Configruation
% for RRR = 1:8 % Number of cties
%     clearvars -except RRR;
%     store = open('file_name.mat');
%     warning off
%     open_name = store.file_name{RRR};
%     store = open(open_name);
% 
%     V = store.V;
%     N = store.N;
%     b = store.b;
%     rho = store.rho;
%     x_nominal = store.x_nominal;
%     y_nominal = store.y_nominal';
%     x_robust = store.store.x;
%     y_robust = store.y_robust;
% 
%     for v = 1:V
%         recourse_nominal(v) = sum(b.*y_nominal{v}, 'all');    
%         recourse_robust(v) = sum(b.*y_robust{v}, 'all');
%     end
%         
%     benefit_abs = recourse_robust - recourse_nominal;
%     benefit_r = 1000*benefit_abs./recourse_nominal;
%     benefit_r = [benefit_r' (1:V)'];
%     benefit_r = sortrows(benefit_r, 'descend');
%     sample_index = benefit_r(:, 2);
%     sample_index = sample_index(1:N);
%     benefit_r = benefit_r(:, 1);
%     benefit_r = benefit_r(1:N);
% 
%     max_r = max(benefit_r);
%     mean_r = mean(benefit_r);
% 
%     benefit_abs = benefit_abs(sample_index);
%     max_abs = max(benefit_abs);
%     mean_abs = mean(benefit_abs);
% 
%     filename = ['Benefit_N' num2str(N) '_Rho' num2str(rho) '.mat'];
%     save(filename);
% 
% end

rho = 0.8;

for RRR = 1:8
    
    b = b_all{RRR};
    x_robust = x_robust_all{RRR};
    y_robust = y_robust_all{RRR};
    x_nominal = x_nominal_all{RRR};
    y_nominal = y_nominal_all{RRR};
    V = V_all(RRR);
    N = N_all(RRR);

    for v = 1:V
        recourse_nominal(v) = sum(b.*y_nominal{v}, 'all');    
        recourse_robust(v) = sum(b.*y_robust{v}, 'all');
    end
     
    benefit_abs = recourse_robust - recourse_nominal;
    benefit_r = 1000*benefit_abs./recourse_nominal;
    benefit_r = [benefit_r' (1:V)'];
    benefit_r = sortrows(benefit_r, 'descend');
    sample_index = benefit_r(:, 2);
    sample_index = sample_index(1:N);
    benefit_r = benefit_r(:, 1);
    benefit_r = benefit_r(1:N);

    max_r = max(benefit_r);
    mean_r = mean(benefit_r);

    benefit_abs = benefit_abs(sample_index);
    max_abs = max(benefit_abs);
    mean_abs = mean(benefit_abs);

    filename = ['Benefit_N' num2str(N) '_Rho' num2str(rho) '.mat'];
    save(filename);

end