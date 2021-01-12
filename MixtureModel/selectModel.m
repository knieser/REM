function [AIC_EM, BIC_EM, AIC_REM, BIC_REM] = selectModel(sim_num,delta)

% Choose number of clusters based on AIC/BIC using REM estimates
    
% Set seed
rng(2021);

load(['GMM_output_Sim_',num2str(sim_num),'_k_2_delta_',num2str(100*delta),'.mat'],'X','global_iter','step'); 
p = size(X,1);
n = size(X,2);

% AIC/BIC EM
AIC_EM = 1:9;
BIC_EM = 1:9;

for j = 1:length(AIC_EM)
    gmmfit = fitgmdist(X',j, 'Options',statset('MaxIter',1500));
    AIC_EM(j) = gmmfit.AIC;
    BIC_EM(j) = gmmfit.BIC;
end

% AIC/BIC REM
AIC_REM = 1:5;
BIC_REM = 1:5;

for k = 1:5

    disp(['Working on k = ', num2str(k)])
    
    [~,~,~,~,~,~,lik] = GMM_estimates_binary_search(X,k,delta,global_iter,step);
    
    parms = k*p + k*p*(p+1)/2 + k-1;
    
    AIC_REM(k) = -2*lik + 2*parms;
    BIC_REM(k) = -2*lik + log(n)*parms;

end

end