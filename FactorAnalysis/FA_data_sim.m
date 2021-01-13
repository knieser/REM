function [sigma_01, sigma_02, lambda_01, lambda_02, R_lambdas] = FA_data_sim(p,k,sparse,communality,m,seed)
%{
This is the function generates the data.
    
INPUT:
    p: number of observed variables
    k: number of latent factors
    sparse: 0/1 indicator
    communality: 1 = low, 2 = wide, 3 = high
    m: values between 0 and 1, degree of similarity between two population
    covariance matrices
    
OUTPUT:
    sigma_01: (p x p) population covariance matrix
    sigma_02: (p x p) population covariance matrix
    lambda_01: (p x k) population loading matrix
    lambda_02: (p x k) population loading matrix
    R_lambdas: congruence between lambda_01 and lambda_02
%}

% Set seed;
rng(seed);

% Generate two different sets of lambdas
% Either sparse structure or complex structure
if sparse == 1
    [sigma_01, lambda_01] = SigmaSim(p,k,communality);
    [sigma_02, lambda_02] = SigmaSim(p,k,communality); 
else
    [sigma_01, lambda_01] = SigmaSimNormals(p,k,communality);
    [sigma_02, lambda_02] = SigmaSimNormals(p,k,communality); 
end

% Control similarility of sigmas; 
sigma_02 = m*sigma_01 + (1-m)*sigma_02;

% Compute congruence coeff between lambda_01 and lambda_02
R_lambdas = computeCongruence(lambda_01, lambda_02);
disp(R_lambdas)    
end