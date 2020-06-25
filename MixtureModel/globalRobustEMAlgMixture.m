function [est_mu, est_sigma, est_mix, est_omega, est_gamma, est_weights, final_llh, final_ind_llh] ...
    = globalRobustEMAlgMixture(X,k,iter,epsilon,robust)
%{
This function runs global optimization for EM and REM estimation.
    
INPUT:
    X: (p x n) dataset
    k: number of latent groups
    iter: number of starting points for global optimization
    epsilon: hyperparameter for REM estimation 
    robust: 0/1 flag for REM estimation
    
OUTPUT:
    est_mu: (p x k) estimated mean vectors
    est_sigma: (p x p x k) estimated covariance matrices
    est_mix: (1 x p) estimated mixture proportions
    est_omega: (n x 1) estimated posterior probabilities
    est_gamma: estimated gamma from REM estimation
    est_weights: individual-level probabilistic weights from REM estimation
    final_llh: joint loglikelihood
    final_ind_llh: (n x 1) individual-level loglikelihood values
%}  
    
% Get data dimensions
p = size(X,1); 
n = size(X,2);
  
% Use MATLAB function to get starting estimates for sigma and mix
GMModel = fitgmdist(X',k,'Options',statset('MaxIter',500));
intl_sigma = GMModel.Sigma;
intl_mix = GMModel.ComponentProportion;

% Random starting points for mu
X_trim = X(:,1:k*floor(n/k));
mu_start = datasample(reshape(X_trim,p,k,[]),iter,3,'Replace',false);

% Initialize likelihood values
max_alt_llh = -Inf;
max_llh = -Inf;
    
if robust == 1   
    for l = 1:iter 
    
        intl_mu = mu_start(:,:,l);
    
        [mu, sigma, mix, omega, gamma, weights, alt_llh, alt_ind_llh] ...
            = RobustEMAlgMixture(X,k,epsilon,intl_mu,intl_sigma,intl_mix);
                
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

        [mu, sigma, mix, omega, llh, ind_llh] = EMAlgMixture(X,k,intl_mu,intl_sigma,intl_mix);
    
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