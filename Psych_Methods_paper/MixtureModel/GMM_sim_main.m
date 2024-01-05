function [X,grp_flag,true_gmm,est_gmm,est_gmm_2,weights,est_gamma,opt_eps,alt_lik,lik] = GMM_sim_main(sim_num,k,delta)
%{
This is the main function for running the simulations of Gaussian mixture
models and obtaining parameter estimates through EM and REM.
    
INPUT:
    sim_num: 
        1 = No corruption,
        2 = Scattered minority,
        3 = Scattered minority w cross
        4 = Scattered minority, greater within-group variance
        5 = Scattered minority, greater within-group variance
        6 = Noise centered within a group
    k: number of latent groups
    delta: hyperparameter for REM estimation 
    
OUTPUT:
    true_gmm: ground truth distribution
    est_gmm: EM estimated distribution
    est_gmm_2: REM estimated distribution
    weights: individual-level probabilistic weights from REM estimation
    est_gamma: estimated gamma from REM estimation
    opt_eps: optimal epsilon hyperparameter for REM estimation
    alt_lik: joint alternative log-likelihood
    lik: joint log-likelihood
%}

% Set seed
rng(2021);

% Set parameters
% Number of starting points for global optimization
global_iter = 500;

% Number of iterations to search for epsilon
step = 100;

% Simulated data parameters
p_0 = 2; % number of dimensions of data
k_0 = 3; % true number of latent groups
n = 1000;
mix = [0.70 0.20 0.10];

% Make data
[X, true_gmm, grp_flag] = GMM_data_sim(sim_num,p_0,k_0,n,mix);

% Get EM and REM estimates
[est_gmm, est_gmm_2, est_weights, est_gamma, opt_eps,alt_lik,lik] = GMM_estimates_binary_search(X,k,delta,global_iter,step);

% Add flag for minority group for weights
minority = [grp_flag grp_flag > k];
weights = [est_weights minority(:,2)];

% Restrict ground truth distribution to estimated latent groups
true_gmm = gmdistribution(true_gmm.mu(1:k,:), true_gmm.Sigma(:,:,1:k), true_gmm.ComponentProportion(1:k));

% Save
output = ['GMM_output_Sim_', num2str(sim_num),'_k_',num2str(k),'_delta_', num2str(100*delta)];
save(output);

end

%%% Calculate sample estimates %%%
% Mean
%   mean(X(:,grp_flag==1),2)
%   mean(X(:,grp_flag==2),2)

% Sigma
%   cov(X(:,grp_flag==1)')
%   cov(X(:,grp_flag==2)')

% Pi
%   sum(grp_flag==1)/(n - sum(grp_flag==3)) 
%   sum(grp_flag==2)/(n - sum(grp_flag==3)) 

