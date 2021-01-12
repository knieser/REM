function [mu, sigma, mix, omega, gamma, weights, alt_llh, ind_llh] ...
        = RobustEMAlgMixture(X,k,epsilon,intl_mu,intl_sigma,intl_mix)
%{
This function obtains REM estimates.
    
INPUT:
    X: (p x n) dataset
    k: number of latent groups
    epsilon: hyperparameter
    intl_mu: (p x k) intial values for mean vectors
    intl_sigma: (p x p x k) initial values for covariance matrices
    intl_mix: (1 x p) initial values for mixture proportions
    
OUTPUT:
    mu: (p x k) estimated mean vectors
    sigma: (p x p x k) estimated covariance matrices
    mix: (1 x p) estimated mixture proportions
    omega: (n x 1) estimated posterior probabilities
    gamma: estimated gamma from REM estimation
    weights: individual-level probabilistic weights from REM estimation
    alt_llh: joint alternative log-likelihood
    ind_llh: (n x 1) individual log-likelihood values
%}  

% Get dimensions
p = size(X,1);
n = size(X,2);

% Initial guesses;
gamma = 0.5;
mu = intl_mu;
sigma = intl_sigma;
mix = intl_mix;
omega = zeros(n,k);

% Tolerance parameters
tol = 1e-15;
maxiter = 500;
llh = -inf(1,maxiter);

% Initialize
log_g = zeros(n,k);
exit = 0; 

% Estimation 
for iter = 2:maxiter 
        
    % Calculate log_g and omega_denom for each cluster;
    omega_denom = zeros(n,1);
    for j = 1:k
        try 
            chol(sigma(:,:,j));      
        catch 
            %warning(['Set gamma to 0 for log(epsilon) = ',num2str(log(epsilon))])
            gamma = 0;
            weights = zeros(n,1);
            alt_ind_llh = -Inf;
            exit = 1;
            break;
        end
        
        V = chol(sigma(:,:,j));
        
        log_g(:,j) = -1/2*(p*log(2*pi) + log(det(sigma(:,:,j))) + diag( ( ((X - mu(:,j))'/V) / (V') ) * (X - mu(:,j)))); 
        % returns n x 1;
        omega_denom = omega_denom + mix(j)*exp(log_g(:,j));
    end
    
    if exit == 1
        break;
    end
    
    % Iteration stop criteria
    ind_llh = log(sum(mix.*exp(log_g),2)); % log-likelihood for each data point;
    alt_ind_llh = log(gamma*sum(mix.*exp(log_g),2) + (1-gamma)*epsilon);
    llh(iter) = sum(alt_ind_llh); % Joint log-likelihood for X;
    if abs(llh(iter) - llh(iter-1)) < tol*abs(llh(iter-1)); break; end
    
    % Gather weights p(y) for robust weighting;
    scaledllh = log(gamma)+ind_llh; % Add in gamma
    weights   = log(1-gamma) + log(epsilon) - scaledllh;
    weights   = 1./( 1 + exp(weights) );

    % Estimate gamma
    gamma = 1/n*sum(weights);
    
    % E Step
    omega = mix.*exp(log_g) ./ omega_denom;     
             
    % M Step
    for j = 1:k
        total_weight = weights.*omega(:,j);
        mu(:,j) = (X * total_weight) ./ sum(total_weight);
        sigma(:,:,j) = (X - mu(:,j))*(total_weight.*(X - mu(:,j))') ./ sum(total_weight);
        mix(:,j) = sum(total_weight) / sum(weights);   
    end 
       
end

alt_llh = sum(alt_ind_llh);

end