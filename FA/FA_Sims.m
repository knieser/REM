function [R_EM, R_REM, gamma, opt_eps, R_EM_minor, R_REM_minor] = FA_Sims(sim)
    
seed = rng();


%%%%% Set Parameters %%%%%%

% mixture rates;
corrupt_pct = 0:0.10:0.40;

% k_choices;
k_choice = 1:3;

% # simulations to run to generate error bars;
n_sim = 10; 

% factor structure parameters;
sparse = 1; % if 1 will use Tucker-Koopman-Linn method to generate pop correlation matrices;
communality = 3; % 1 = low, 2 = wide, 3 = high, 4 = perfect;
m = 0; % controls similarity of lambdas, m=1 means the lambdas are the same;

% dimensions
p_0 = 10;
k_0 = 1;
n = 1000;

% estimation parameters;
delta = 0.02;
step = 100;


%%%%%% Run Simulation %%%%%%

if sim == 1
    msg = 'Simulation 1: Varying Mix';
    disp(msg)
    
    [R_EM, R_REM, gamma, opt_eps, R_EM_minor, R_REM_minor] = VaryMix(p_0,k_0,n,k_choice,corrupt_pct,sparse,communality,m,delta,step,n_sim,seed);

elseif sim == 2
    msg = 'Simulation 2: Varying K';
    disp(msg)
    
    [R_EM, R_REM, gamma, opt_eps, R_EM_minor, R_REM_minor] = VaryK(p_0,k_0,n,k_choice,corrupt_pct,sparse,communality,m,delta,step,n_sim,seed);

end



end



