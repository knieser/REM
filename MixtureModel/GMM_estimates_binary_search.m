function [est_gmm,est_gmm_2,est_weights,opt_gamma,opt_eps,alt_lik,lik] = GMM_estimates_binary_search(X,k,delta,global_iter,step)
%{
This is the main function for obtaining EM and REM estimates.
    
INPUT:
    X: (p x n) dataset
    k: number of latent groups
    delta: hyperparameter for REM estimation 
    global_iter: number of starting points for global optimization
    step: number of iterations to search for epsilon
    
OUTPUT:
    est_gmm: EM estimated distribution
    est_gmm_2: REM estimated distribution
    est_weights: individual-level probabilistic weights from REM estimation
    est_gamma: estimated gamma from REM estimation
    opt_eps: optimal epsilon hyperparameter for REM estimation
    alt_lik: joint alternative log-likelihood
    lik: joint log-likelihood
%}  

%%%% Get estimates from EM %%%%
disp('Starting EM global opt...')

[mu_1, sigma_1, mix_1, ~, ~, ~, ~, ind_llh] ...
    = globalRobustEMAlgMixture(X,k,global_iter,0,0); 

disp('...done')

%%%% Get estimates from REM %%%%

% Get initial seeds
GMModel = fitgmdist(X',k,'Options',statset('MaxIter',500));
intl_sigma = GMModel.Sigma;
intl_mix = GMModel.ComponentProportion;

% Global optimization to get starting mu
eps_start = quantile(exp(ind_llh),0.05);
disp('Searching for REM global opt...')

[intl_mu] = globalRobustEMAlgMixture(X,k,global_iter,eps_start,1);

disp('...done')

% Set max of epsilon to check
eps_max = quantile(exp(ind_llh), 0.2); % Median likelihood value

% Initialize lower and upper bounds of search for epsilon
eps_range = zeros(step,1);
eps_range(2) = eps_max;
est_gamma = zeros(length(eps_range),1);
chk = zeros(length(eps_range),1);

% Check lower bound
leps = 0;
est_gamma(1) = 1;
chk(1) = 0;

% Check upper bound
ueps = eps_max;
[mu_2, sigma_2, mix_2, ~, est_gamma(2)] = RobustEMAlgMixture(X,k,ueps,intl_mu,intl_sigma,intl_mix);
chk(2) = checkEps(mu_2,sigma_2,mix_2,0.90,eps_range(2));

% Update upper bound if needed
while chk(2) < delta
    
    % If gamma is too small, stop increasing upper bound
    if est_gamma(2) < 0.5
        break;
    end
    
    % Display
    disp('Expanding search range for epsilon')
        
    % Extend upper bound
    eps_range(2) = eps_range(2)*2;
    ueps = eps_range(2);
    [mu_2, sigma_2, mix_2, ~, est_gamma(2)] = RobustEMAlgMixture(X,k,ueps,intl_mu,intl_sigma,intl_mix);
    chk(2) = checkEps(mu_2,sigma_2,mix_2,0.90,eps_range(2));
    
end

% Set next epsilon to check
eps_range(3) = (leps+ueps)/2;

% Run through binary search of epsilon
for iter = 3:step

    % Estimate parameters with epsilon
    [mu_2, sigma_2, mix_2, ~, est_gamma(iter)] = RobustEMAlgMixture(X,k,eps_range(iter),intl_mu,intl_sigma,intl_mix);
    
    % Check heuristic for hyperparameter selection;
    chk(iter) = checkEps(mu_2,sigma_2,mix_2,0.90,eps_range(iter)); %est_gamma(iter)
            
    % Update upper and lower bounds
    if chk(iter) < delta
        % Increase lower bound (account for noise)
        leps = leps*(1/10) + eps_range(iter)*(9/10);
    else
        % Decrease upper bound (account for noise)
        ueps = ueps*(1/10) + eps_range(iter)*(9/10);
    end
    
    % Set new epsilon to check
    eps_range(iter+1) = (leps+ueps)/2;
    
end

% Store values;
est_gmm = gmdistribution(mu_1', sigma_1, mix_1);
opt_eps = eps_range(iter+1);
[mu_2, sigma_2, mix_2, ~, opt_gamma,est_weights,alt_lik,ind_llh] = RobustEMAlgMixture(X,k,opt_eps,intl_mu,intl_sigma,intl_mix);
lik = sum(ind_llh);
est_gmm_2 = gmdistribution(mu_2', sigma_2, mix_2);

end