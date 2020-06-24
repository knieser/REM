function [true_gmm,est_gmm,est_gmm_2,weights,est_gamma,opt_eps] = GMM_sim_main(sim_num,k,delta)
% Sim Key
% 1 = No corruption
% 2 = Scattered minority
% 3 = Scattered minority w cross

% e.g. LPA_sim_main(1,2,0.05);

seed = rng(444);

%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%

% Number of starting points for global optimization;
global_iter = 500;

% Number of iterations to search for epsilon;
step = 100;

% Simulated data parameters;
p_0 = 2; % number of dimensions of data;
k_0 = 3; % true number of clusters;
n = 1000;
mix = [0.70 0.20 0.10];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make data
[X, true_gmm, grp_flag] = GMM_data_sim(sim_num,p_0,k_0,n,mix,seed);

% Compute AIC and BIC
AIC = 1:k_0+1;
BIC = 1:k_0+1;

for j = 1:5
    gmmfit = fitgmdist(X',j);
    AIC(j) = gmmfit.AIC;
    BIC(j) = gmmfit.BIC;
end

% Get estimates
[est_gmm, est_gmm_2, est_weights, est_gamma, opt_eps] = GMM_estimates(X,k,delta,global_iter,step,seed);

% Add flag for minority group for weights;
minority = [grp_flag grp_flag > k];
weights = [est_weights minority(:,2)];

% Limit to clusters we are interested in for figure
true_gmm = gmdistribution(true_gmm.mu(1:k,:), true_gmm.Sigma(:,:,1:k), true_gmm.ComponentProportion(1:k));

% Save
output = ['GMM_output_Sim_', num2str(sim_num),'_k_',num2str(k),'_delta_', num2str(100*delta)];
save(output);

end

% Calculate sample estimates;
% Mean
% mean(X(:,grp_flag==1),2)
% mean(X(:,grp_flag==2),2)

% Sigma
% cov(X(:,grp_flag==1)')
% cov(X(:,grp_flag==2)')

% Pi
% sum(grp_flag==1)/(n - sum(grp_flag==3)) 
% sum(grp_flag==2)/(n - sum(grp_flag==3)) 

