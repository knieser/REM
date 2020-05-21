function [true_gmm,est_gmm,est_gmm_2,weights,est_gamma,opt_eps] = LPA_sim_main(sim_num,k,delta)
% Sim Key
% 1 = No corruption
% 2 = Scattered minority
% 3 = Scattered minority w cross

% e.g. LPA_sim_main(1,2,0.05);
  
%%%%%%%%%%% Parameters %%%%%%%%%%%%%%%

% Number of starting points for global optimization;
global_iter = 300;

% Number of iterations to search for epsilon;
step = 100;

% Simulated data parameters;
p_0 = 2; % number of dimensions of data;
k_0 = 3; % true number of clusters;
n = 1000;
mix = [0.70 0.20 0.10];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make data
[X, true_gmm, grp_flag] = LPA_data_sim(sim_num,p_0,k_0,n,mix);

% Compute AIC and BIC
AIC = 1:k_0+1;
BIC = 1:k_0+1;

for j = 1:k_0+1
    gmmfit = fitgmdist(X',j);
    AIC(j) = gmmfit.AIC;
    BIC(j) = gmmfit.BIC;
end

% Get estimates
[est_gmm, est_gmm_2, est_weights, est_gamma, opt_eps] = LPA_estimates(X,k,delta,global_iter,step);

% Add flag for minority group for weights;
minority = [grp_flag grp_flag > k];
weights = [est_weights minority(:,2)];

% Limit to clusters we are interested in for figure
true_gmm = gmdistribution(true_gmm.mu(1:k,:), true_gmm.Sigma(:,:,1:k), true_gmm.ComponentProportion(1:k));

% Save
output = ['LPA_output_Sim_', num2str(sim_num),'_k_',num2str(k),'_delta_', num2str(100*delta)];
save(output);

end