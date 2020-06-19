function [sigma_01, sigma_02, lambda_01, lambda_02, R_sigmas] = FA_data_sim(p,k,sparse,communality,m,seed)

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

% Control similarility of lambdas; m=1 means the correlation matrices are the same
sigma_02 = m*sigma_01 + (1-m)*sigma_02;

% Compute congruence coeff between lambda_01 and lambda_02
R_sigmas = computeCongruence(lambda_01, lambda_02);
    
end