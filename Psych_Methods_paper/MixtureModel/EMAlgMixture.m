function [mu, sigma, mix, omega, llh, ind_llh] = EMAlgMixture(X,k,mu,sigma,mix)
%{
This function obtains EM estimates.
    
INPUT:
    X: (p x n) dataset
    k: number of latent groups
    mu: (p x k) intial values for mean vectors
    sigma: (p x p x k) initial values for covariance matrices
    mix: (1 x p) initial values for mixture proportions
    
OUTPUT:
    mu: (p x k) estimated mean vectors
    sigma: (p x p x k) estimated covariance matrices
    mix: (1 x p) estimated mixture proportions
    omega: (n x 1) estimated posterior probabilities
    llh: joint log-likelihood
    ind_llh: (n x 1) individual log-likelihood values
%}  

% Get dimensions
p = size(X,1);
n = size(X,2);

% Tolerance parameters
tol = 1e-15;
maxiter = 500;
llh = -inf(1,maxiter);

% Initialize
log_g = zeros(n,k);

% Estimation 
for iter = 2:maxiter
    
    % Calculate log_g and omega_denom for each cluster;
    omega_denom = zeros(n,1);
    for j = 1:k
        
        try 
            chol(sigma(:,:,j));      
        catch 
            llh = -Inf;
            return;
        end
        
        log_g(:,j) = -1/2*(p*log(2*pi) + log(det(sigma(:,:,j))) + diag((X - mu(:,j))'/sigma(:,:,j)*(X - mu(:,j)))); 
        % returns n x 1;
        omega_denom = omega_denom + mix(j)*exp(log_g(:,j));
    end
    
    % Iteration stop criteria
    ind_llh = log(sum(mix.*exp(log_g),2));
    llh(iter) = sum(ind_llh); % Joint loglikelihood for X;
    if abs(llh(iter) - llh(iter-1)) < tol*abs(llh(iter-1)); break; end
    
    % E Step
    omega = mix.*exp(log_g) ./ omega_denom;     
             
    % M Step
    for j = 1:k
        mu(:,j) = (X * omega(:,j)) ./ sum(omega(:,j));
        sigma(:,:,j) = (X - mu(:,j))*(omega(:,j).*(X - mu(:,j))') ./ sum(omega(:,j));
        mix(:,j) = 1/n*sum(omega(:,j));   
    end 
       
end

llh = llh(iter);

end