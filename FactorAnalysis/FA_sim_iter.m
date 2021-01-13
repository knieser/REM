function [R_EM, R_REM, R_EM_minor, R_REM_minor, gamma_values, eps_values, err_values] = ...
    FA_sim_iter(n_sim,sigma_01,sigma_02,lambda_01,lambda_02,n,corrupt_pct,k,delta)
%{
This is the function that samples data and obtains aggregated EM and REM
estimates for a given level of corruption and number of latent factors
    
INPUT:
    n_sim: number of simulations to run
    sigma_01: (p x p) population covariance matrix
    sigma_02: (p x p) population covariance matrix
    lambda_01: (p x k) population loading matrix
    lambda_02: (p x k) population loading matrix
    n: number of observations
    corrupt_pct: proportion of minority group  
    k: number of latent factors
    delta: hyperparameter for REM estimation
    
OUTPUT:
    R_EM: mean and SE of congruence coeff of EM estimate with majority factor
    structure
    R_REM: mean and SE of congruence coeff of REM estimate with majority factor
    structure
    R_EM_minor: mean and SE of congruence coeff of EM estimate with minority factor
    structure
    R_REM_minor: mean and SE of congruence coeff of REM estimate with minority factor
    structure
    gamma_values: mean and SE of gamma
    eps_values: mean and SE of epsilon
%}

% Initialize
R_classic = zeros(n_sim,1);
R_robust = zeros(n_sim,1);
R_classic_b = zeros(n_sim,1);
R_robust_b = zeros(n_sim,1);
gamma = zeros(n_sim,1);
eps = zeros(n_sim,1);
err = zeros(n_sim,1);

parfor l = 1:n_sim   
    % Display iteration number
    if mod(l-1,10)==0 
        disp(['Iteration ',num2str(l),' of ',num2str(n_sim)])
    end
    
    % Sample data from population;
    try
    [X,lambda_target_1,lambda_target_2] ...
        = FA_sample_data(sigma_01, sigma_02, lambda_01, lambda_02,n,corrupt_pct,l);
    
    % Get EM and REM estimates;
    [hlambda_1, ~, hlambda_2, ~,  ~, gamma(l), eps(l), err(l)] = FA_estimates(X,k,delta);
    
    % Compute and store congruence coefficient
    R_classic(l) = computeCongruence(hlambda_1, lambda_target_1); % first matrix gets procrustes rotation;
    R_robust(l) = computeCongruence(hlambda_2, lambda_target_1);
    
    R_classic_b(l) = computeCongruence(hlambda_1, lambda_target_2); % first matrix gets procrustes rotation;
    R_robust_b(l) = computeCongruence(hlambda_2, lambda_target_2);   
    catch 
       disp('ERROR!! Skipping iteration') 
    end
end

% Aggregate data;
R_EM = aggregate(R_classic);
R_REM = aggregate(R_robust);

R_EM_minor = aggregate(R_classic_b);
R_REM_minor = aggregate(R_robust_b);

gamma_values = aggregate(gamma);
eps_values = aggregate(eps);
err_values = aggregate(err);

end