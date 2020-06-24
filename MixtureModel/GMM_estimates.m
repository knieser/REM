function [est_gmm,est_gmm_2,est_weights,est_gamma,opt_eps] = GMM_estimates(X,k,delta,global_iter,step,seed)
    
rng(seed);

% Set # iterations for epsilon search;
maxiter = step+1;
    
% Get estimates from EM; 
disp('Starting EM global opt...')
[mu_1, sigma_1, mix_1, ~, ~, ~, ~, ind_llh] ...
    = globalRobustEMAlgMixture(X,k,global_iter,0,0,seed);  
disp('...done')

% Get estimates from REM;
chk = zeros(1,maxiter);
est_gamma = zeros(1,maxiter);

eps_max = quantile(exp(ind_llh), 0.50); % Median likelihood value;
eps_range = 0:eps_max/step:eps_max;

% Get initial seeds;
GMModel = fitgmdist(X',k,'Options',statset('MaxIter',500));
intl_sigma = GMModel.Sigma;
intl_mix = GMModel.ComponentProportion;

% Global optimization to get starting mu;
eps_start = quantile(exp(ind_llh),0.05);
disp('Searching for REM global opt...')
[intl_mu] = globalRobustEMAlgMixture(X, k, global_iter,eps_start,1,seed);
disp('...done')

for iter = 1:maxiter

    % Get REM estimates;
    [mu_2, sigma_2, mix_2, ~, est_gamma(iter), est_weights] ...
        = RobustEMAlgMixture(X,k,eps_range(iter),intl_mu,intl_sigma,intl_mix);
    
    % Calculate prob. that weights <= 1/2;
    chk(iter) = checkEps(mu_2,sigma_2,mix_2,est_gamma(iter),eps_range(iter),seed);
            
    if chk(iter) > delta 
        break;
    end
end
        
opt_eps = eps_range(iter);
est_gmm = gmdistribution(mu_1', sigma_1, mix_1);
est_gmm_2 = gmdistribution(mu_2', sigma_2, mix_2);

end