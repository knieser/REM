function [X,lambda_target_1,lambda_target_2, grp_flag] = ...
    FA_sample_data(sigma_01, sigma_02, lambda_01, lambda_02,n,corrupt_pct)
%{
This is the function samples data from population.
    
INPUT:
    sigma_01: (p x p) population covariance matrix
    sigma_02: (p x p) population covariance matrix
    lambda_01: (p x k) population loading matrix
    lambda_02: (p x k) population loading matrix
    n: number of observations
    corrupt_pct: proportion of minority group  
    
OUTPUT:
    X: (p x n) array of simulated data
    lambda_target_1: normalized population loading matrix for majority
    lambda_target_2: normalized population loading matrix for minority
    grp_flag: nx1 vector of group indicators
%}
    
p = size(sigma_01,1);

corrupt_number = ceil(n*(corrupt_pct));

% Sample data from each population
R_01 = chol(sigma_01);
X_01 = randn(n,p)*R_01;
X_01 = X_01';

R_02 = chol(sigma_02);
X_02 = randn(n,p)*R_02;
X_02 = X_02';

% Create flag majority and minority groups
grp_flag = [ones(1,n-corrupt_number) zeros(1,corrupt_number)];

% Mix samples
X = [X_01(:,1:n-corrupt_number) X_02(:,1:corrupt_number)];  

% Normalize data and lambdas;
V = var(X,0,2);
invV = 1./V;
X = normalize(X')';
lambda_target_1 = sqrt(invV).*lambda_01;
lambda_target_2 = sqrt(invV).*lambda_02;
    
end