function [hlambda_1, hpsi_1, hlambda_2, hpsi_2, weights, gamma, eps] = FA_estimates(X,k,delta)
%{
This is the main function for obtaining EM and REM estimates.
    
INPUT:
    X: (p x n) dataset
    k: number of latent factors
    delta: hyperparameter for REM estimation 
    
OUTPUT:
    hlambda_1: (p x k) EM estimated lambda
    hpsi_1: (p x 1) EM estimated psi
    hlambda_2: (p x k) REM estimated lambda
    hpsi_2: (p x 1) REM estimated psi
    weights: individual-level probabilistic weights from REM estimation
    gamma: estimated gamma from REM estimation
    eps: optimal epsilon hyperparameter for REM estimation
%}  
    
% Number of items    
p = size(X,1);

% Set # iterations for epsilon search
step = 100;

% Get EM estimates
[hlambda_1, hpsi_1, ind_llh] = EMAlg(X,k);    
            
% Get REM estimates
intl_lambda = hlambda_1 + randn(p,k);
intl_psi = hpsi_1;

% Set range of epsilon to check;
eps_max = quantile(exp(ind_llh), 0.50); % Median likelihood value;
eps_range = 0:eps_max/step:eps_max;

est_gamma = zeros(length(eps_range));
chk = zeros(length(eps_range));

for iter = 1:length(eps_range)
    
    [hlambda_2, hpsi_2, est_gamma(iter), weights] = RobustEMAlg(X,k,eps_range(iter),intl_lambda,intl_psi);
                    
    % Check heuristic for hyperparameter selection
    chk(iter) = checkEpsFA(hlambda_2,hpsi_2,est_gamma(iter),eps_range(iter));
    
    if chk(iter) > delta 
        break;
    end
end
            
% Store values
gamma = est_gamma(iter);  
eps = eps_range(iter);
     
end

