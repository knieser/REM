function [R_EM, R_REM, gamma_values, eps_values, R_EM_minor, R_REM_minor] = FA_sim_main(sim_num,communality,delta)

% sim: (1) vary corruption (2) vary K
% can enter array into either corrupt_pct or k_choice
%communality: 1 = low, 2 = wide, 3 = high;

% e.g. FA_sim_main(1,0:0.1:0.4,1,3,0.05);

seed = rng(123);

%%%%%%%%% Parameters %%%%%%%%%%%%

% # simulations to run to generate error bars;
n_sim = 100; 

% factor structure parameters;
sparse = 1; % if 1 will use Tucker-Koopman-Linn method to generate pop correlation matrices;
m = 0; % controls similarity of lambdas, m=1 means the lambdas are the same;

% dimensions
n = 600;
if sim_num == 1
    p_0 = 30;
    k_0 = 4;
    corrupt_pct = 0:0.05:0.4;
    k_choice = k_0;
elseif sim_num == 2
    p_0 = 30;
    k_0 = 4;
    corrupt_pct = 0.30;
    k_choice = 1:k_0+2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Simulate population covariance/factor structures 
[sigma_01, sigma_02, lambda_01, lambda_02, ~] = FA_data_sim(p_0,k_0,sparse,communality,m,seed);

% Calculate communalities
pop_h2 = sum(lambda_01.^2,2); % p x 1;

% Run simulation 
if sim_num == 1
    msg = 'Simulation 1: Varying Mix';
    disp(msg)
    
    % Initialize;
    R_EM = zeros(length(corrupt_pct), 2);
    R_REM = zeros(length(corrupt_pct), 2);
    
    R_EM_minor = zeros(length(corrupt_pct), 2);
    R_REM_minor = zeros(length(corrupt_pct), 2);
    
    gamma_values = zeros(length(corrupt_pct), 2);
    eps_values = zeros(length(corrupt_pct), 2);
    
    avg_EM_h2 = zeros(p_0,length(corrupt_pct));
    avg_REM_h2 = zeros(p_0,length(corrupt_pct));
    
    % Run through varying corrupt_pct
    for c = 1:length(corrupt_pct)
        msg = ['Working on corrupt pct = ', num2str(100*corrupt_pct(c)),'%'];
        disp(msg)
        
        [R_EM(c,:), R_REM(c,:), R_EM_minor(c,:), R_REM_minor(c,:), avg_EM_h2(:,c), avg_REM_h2(:,c), gamma_values(c,:), eps_values(c,:)] = ... 
        FA_sim_iter(n_sim,sigma_01,sigma_02,lambda_01,lambda_02,n,corrupt_pct(c),k_choice,delta,seed);
    end 

    output = ['FA_output_Sim_', num2str(sim_num),'_comm_',num2str(communality),'_delta_', num2str(100*delta),'_test'];
    
elseif sim_num == 2
    msg = 'Simulation 2: Varying K';
    disp(msg)
    
    % Initialize;
    R_EM = zeros(length(k_choice), 2);
    R_REM = zeros(length(k_choice), 2);
    
    R_EM_minor = zeros(length(k_choice), 2);
    R_REM_minor = zeros(length(k_choice), 2);
    
    gamma_values = zeros(length(k_choice), 2);
    eps_values = zeros(length(k_choice), 2);
    
    avg_EM_h2 = zeros(p_0,length(k_choice));
    avg_REM_h2 = zeros(p_0,length(k_choice));
         
    % Run through varying k_choice
    for k = 1:length(k_choice)
        msg = ['Working on k = ', num2str(k_choice(k))];
        disp(msg)
        
        [R_EM(k,:), R_REM(k,:), R_EM_minor(k,:), R_REM_minor(k,:), avg_EM_h2(:,k), avg_REM_h2(:,k), gamma_values(k,:), eps_values(k,:)] = ... 
        FA_sim_iter(n_sim,sigma_01,sigma_02,lambda_01,lambda_02,n,corrupt_pct,k_choice(k),delta,seed);
    end 
    
    output = ['FA_output_Sim_', num2str(sim_num),'_corrupt_',num2str(100*corrupt_pct),'pct_delta_', num2str(100*delta)];
end

save(output);

end



