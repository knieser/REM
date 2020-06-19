function [hlambda_1, hpsi_1, hlambda_2, hpsi_2, weights, gamma, eps] = FA_estimates(X,k,delta)

p = size(X,1);

% Set # iterations for epsilon search;
step = 100;

% Get EM estimates;
[hlambda_1, hpsi_1, ind_llh] = EMAlg(X,k);    
            
% Get REM estimates;
intl_lambda = hlambda_1 + randn(p,k);
intl_psi = hpsi_1;

% Set range of epsilon to check;
eps_max = quantile(exp(ind_llh), 0.50); % Median likelihood value;
eps_range = 0:eps_max/step:eps_max;

est_gamma = zeros(length(eps_range));
chk = zeros(length(eps_range));

for iter = 1:length(eps_range)
    
    [hlambda_2, hpsi_2, est_gamma(iter), weights] = RobustEMAlg(X,k,eps_range(iter),intl_lambda,intl_psi);
                    
    % Calculate prob. that weights <= 1/2;
    chk(iter) = checkEpsFA(hlambda_2,hpsi_2,est_gamma(iter),eps_range(iter));
    
    if chk(iter) > delta 
        break;
    end
end
            
% Store gamma value
gamma = est_gamma(iter);  

% Store epsilon choice
eps = eps_range(iter);
    
    
end

