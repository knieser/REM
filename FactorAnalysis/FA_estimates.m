function [hlambda_1, hpsi_1, hlambda_2, hpsi_2, opt_weights, opt_gamma, opt_eps, opt_error] = FA_estimates(X,k,delta)
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
step = 10;

% Get EM estimates
[hlambda_1, hpsi_1, ind_llh] = EMAlg(X,k);    
            
%%%% Get REM estimates %%%%
intl_lambda = hlambda_1 + randn(p,k);
intl_psi = hpsi_1;

% Set max of epsilon to check;
eps_max   = quantile(exp(ind_llh), 0.22); % Max epsilon value;

% Initialize lower and upper bounds of search for epsilon
eps_range    = zeros(step,1);
eps_range(2) = eps_max;
est_gamma    = zeros(length(eps_range),1);
chk          = zeros(length(eps_range),1);

% Check lower bound
leps         = 0;
est_gamma(1) = 1;
chk(1)       = 0;

% Check upper bound
ueps                              = eps_max;
[hlambda_2, hpsi_2, est_gamma(2)] = RobustEMAlg(X,k,ueps,intl_lambda,intl_psi);
chk(2)                            = checkEpsFA(hlambda_2,hpsi_2,eps_range(2),0.9); 

% Shrink epsilon if already too large
while est_gamma(2) < 0.5
    ueps                              = eps_range(2)/2;
    eps_range(2)                      = ueps;
    [hlambda_2, hpsi_2, est_gamma(2)] = RobustEMAlg(X,k,ueps,intl_lambda,intl_psi);
    chk(2)                            = checkEpsFA(hlambda_2,hpsi_2,eps_range(2),0.9); 
end

% Update upper bound if needed
last_chk = 0;
while chk(2) < delta && last_chk < chk(2) && est_gamma(2) >= 0.5

     % update
     last_chk = chk(2); 
    
     % Extend upper bound
     eps_range(2) = eps_range(2)*2;
     ueps         = eps_range(2);
     [hlambda_2, hpsi_2, est_gamma(2)] = RobustEMAlg(X,k,ueps,intl_lambda,intl_psi);
     chk(2)                            = checkEpsFA(hlambda_2,hpsi_2,eps_range(2),0.9); 
   
end

% Error
if last_chk >= chk(2) ||  est_gamma(2) < 0.5
     disp(['Stopped expanding epsilon range, since stopped increasing'])
     eps_range(2) = eps_range(2)/2;
     ueps         = eps_range(2);
     [hlambda_2, hpsi_2, est_gamma(2)] = RobustEMAlg(X,k,ueps,intl_lambda,intl_psi);
     chk(2)                            = checkEpsFA(hlambda_2,hpsi_2,eps_range(2),0.9);  
end
    
% Set next epsilon to check
eps_range(3) = (leps+ueps)/2;

% Run through binary search of epsilon
for iter = 3:step
    
    % Estimate parameters with epsilon
    [hlambda_2, hpsi_2, est_gamma(iter)] = RobustEMAlg(X,k,eps_range(iter),intl_lambda,intl_psi);
                    
    % Check heuristic for hyperparameter selection
    chk(iter) = checkEpsFA(hlambda_2,hpsi_2,eps_range(iter),0.9); 
    
    % Update upper and lower bounds
    if chk(iter) < delta 
        % Increase lower bound (account for noise)
        leps   = leps*(1/10) + eps_range(iter)*(9/10);
    else
        % Decrease upper bound (account for noise)
        ueps   = ueps*(1/10) + eps_range(iter)*(9/10);
    end
    
    % Set new epsilon to check
    eps_range(iter+1) = (leps+ueps)/2;

end
            
% Store values
opt_eps = eps_range(iter+1);

% Fit REM to full data set using chosen epsilon;
[hlambda_2, hpsi_2, opt_gamma, opt_weights, ~, opt_error] = RobustEMAlg(X,k,opt_eps,intl_lambda,intl_psi);


end

