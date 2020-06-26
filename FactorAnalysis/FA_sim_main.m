function [R_EM, R_REM, gamma_values, eps_values, R_EM_minor, R_REM_minor] ...
        = FA_sim_main(sim_num,communality,delta)
%{
This is the main function for running the simulations of factor structures
and obtaining parameter estimates through EM and REM.
    
INPUT:
    sim_num: 
        1 = Varying mixture percentages,
        2 = Varying number of specified latent factors
    communality: 1 = low, 2 = wide, 3 = high;
    delta: hyperparameter for REM estimation 
    
OUTPUT:
    R_EM: mean and SE of congruence coeff of EM estimate with majority factor
    structure
    R_REM: mean and SE of congruence coeff of REM estimate with majority factor
    structure
    gamma_values: mean and SE of gamma
    eps_values: mean and SE of epsilon
    R_EM_minor: mean and SE of congruence coeff of EM estimate with minority factor
    structure
    R_REM_minor: mean and SE of congruence coeff of REM estimate with minority factor
    structure
%}

% Set parameters

% # simulations to run to generate error bars;
n_sim = 100; 

% Factor structure parameters;
sparse = 1; % if 1 will use Tucker-Koopman-Linn method to generate pop correlation matrices;
m = 0; % controls similarity of lambdas, m=1 means the lambdas are the same;

% Dimensions
n = 600;
p_0 = 30;
k_0 = 4;

if sim_num == 1
    corrupt_pct = 0:0.05:0.4;
    k_choice = k_0;
elseif sim_num == 2
    corrupt_pct = 0.30;
    k_choice = 1:k_0+2;
end

% Simulate population covariance/factor structures 
[sigma_01, sigma_02, lambda_01, lambda_02, ~] = FA_data_sim(p_0,k_0,sparse,communality,m,seed);

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
    
    % Run through varying corrupt_pct
    for c = 1:length(corrupt_pct)
        msg = ['Working on corrupt pct = ', num2str(100*corrupt_pct(c)),'%'];
        disp(msg)
        
        [R_EM(c,:), R_REM(c,:), R_EM_minor(c,:), R_REM_minor(c,:), gamma_values(c,:), eps_values(c,:)] = ... 
        FA_sim_iter(n_sim,sigma_01,sigma_02,lambda_01,lambda_02,n,corrupt_pct(c),k_choice,delta);
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
         
    % Run through varying k_choice
    for k = 1:length(k_choice)
        msg = ['Working on k = ', num2str(k_choice(k))];
        disp(msg)
        
        [R_EM(k,:), R_REM(k,:), R_EM_minor(k,:), R_REM_minor(k,:), gamma_values(k,:), eps_values(k,:)] = ... 
        FA_sim_iter(n_sim,sigma_01,sigma_02,lambda_01,lambda_02,n,corrupt_pct,k_choice(k),delta);
    end 
    
    output = ['FA_output_Sim_', num2str(sim_num),'_corrupt_',num2str(100*corrupt_pct),'pct_delta_', num2str(100*delta)];
end

save(output);

end



