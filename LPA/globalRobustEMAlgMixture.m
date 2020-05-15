function [est_mu, est_sigma, est_mix, est_omega, est_gamma, est_weights, final_llh, final_ind_llh] = globalRobustEMAlgMixture(X,k,iter,epsilon,robust)
  
% Get dimension
p = size(X,1); 
n = size(X,2);
  
% Get data range
%data_rng = [min(X,[],2) max(X,[],2)]; % get range for each dimension
%data_rng = [min(data_rng(:,1)) max(data_rng(:,2))]; % take smallest min and largest max so all data are covered;
 
% Use matlab function to get starting estimates for sigma and mix;
GMModel = fitgmdist(X',k,'Options',statset('MaxIter',500));
intl_sigma = GMModel.Sigma;
intl_mix = GMModel.ComponentProportion;


% Random starting points for mu - uniform over range of data;
%mu_start = data_rng(1) + (data_rng(2) - data_rng(1))*rand(p,k,iter);
X_trim = X(:,1:k*floor(n/k));
mu_start = datasample(reshape(X_trim,p,k,[]),iter,3,'Replace',false);

max_alt_llh = -Inf;
max_llh = -Inf;
    
if robust == 1   
    for l = 1:iter 
    
        intl_mu = mu_start(:,:,l);
    
        [mu, sigma, mix, omega, gamma, weights, alt_llh, alt_ind_llh] = RobustEMAlgMixture(X,k,epsilon, intl_mu, intl_sigma, intl_mix);
                
        if alt_llh > max_alt_llh 
            max_alt_llh = alt_llh;
            est_mu = mu;
            est_sigma = sigma;
            est_mix = mix;
            est_omega = omega;
            est_gamma = gamma;
            est_weights = weights;
            est_ind_llh = alt_ind_llh;
        end
    end 
    final_llh = max_alt_llh;
    final_ind_llh = est_ind_llh;
else
    for l = 1:iter 

        intl_mu = mu_start(:,:,l);

        [mu, sigma, mix, omega, llh, ind_llh] = EMAlgMixture(X, k, intl_mu, intl_sigma, intl_mix);
    
        if llh > max_llh 
            max_llh = llh;
            est_mu = mu;
            est_sigma = sigma;
            est_mix = mix;
            est_omega = omega;
            est_ind_llh = ind_llh;
        end
    end 
    est_gamma = 1;
    est_weights = 1;
    final_llh = max_llh;
    final_ind_llh = est_ind_llh;
end

    
end