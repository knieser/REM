function [R_EM, R_REM, gamma, opt_eps, R_EM_minor, R_REM_minor] = FA_Hetero_Pop_Sim(p_0,k_0,n,k,corrupt_pct,sparse,communality,m,delta,step,n_sim,seed)

rng(seed)
    
% Set # iterations for epsilon search;
maxiter = step+1;

% Choice for p,k in estimation;
p = p_0;

corrupt_number = ceil(n*(corrupt_pct));


%%%%%%%%%%%%% Initialize %%%%%%%%%%%%%%%%%%%%%%%%%%

R_classic = zeros(n_sim,1);
R_classic_b = zeros(n_sim,1);

R_robust = zeros(n_sim,1);
R_robust_b = zeros(n_sim,1);

gamma_values = ones(n_sim,1);
opt_eps = zeros(n_sim,1);



%%%%%%%%%% Data Simulation %%%%%%%%%%%%%%%%%%%%

% Generate two different sets of lambdas
% Either sparse structure or complex structure
if sparse == 1
    [sigma_01, lambda_01] = SigmaSim(p_0,k_0,communality);
    [sigma_02, lambda_02] = SigmaSim(p_0,k_0,communality); 
else
    [sigma_01, lambda_01] = SigmaSimNormals(p_0,k_0,communality);
    [sigma_02, lambda_02] = SigmaSimNormals(p_0,k_0,communality); 
end

% Control similarility of lambdas; m=1 means the correlation matrices are the same
sigma_02 = m*sigma_01 + (1-m)*sigma_02;

% Compute congruence coeff between lambda_01 and lambda_02
R_sigmas = computeCongruence(lambda_01, lambda_02);


%%%%%%%%% Estimation %%%%%%%%%%%%%

for l = 1:n_sim
 
    % Sample data from multivariate normal with mean 0 and covariance sigma;
    R_01 = chol(sigma_01);
    X_01 = randn(n,p)*R_01;
    X_01 = X_01';
    
    R_02 = chol(sigma_02);
    X_02 = randn(n,p)*R_02;
    X_02 = X_02';

    X = [X_01(:,1:n-corrupt_number) X_02(:,1:corrupt_number)];  
    mu = mean(X,2);
    V = var(X,0,2);
    
    X = normalize(X')';

    % Scale lambda_01 to align with normalized data;
    invV = 1./V;
    scaled_lambda_01 = sqrt(invV).*lambda_01;
    scaled_lambda_02 = sqrt(invV).*lambda_02;
    
     % Just focus on the "main" lambda at this point (ie lambda_01);
    lambda_t = scaled_lambda_01;  
    lambda_b = scaled_lambda_02;
        
    %%%%% Get estimates for each sample and mixture pct %%%%%
    
    % Get estimates from EM algorithm;
    [hlambda_1, hpsi_1, ind_llh] = EMAlg(X,k);
    
    chk = zeros(1,maxiter);
    est_gamma = zeros(1,maxiter);
    
    eps_max = quantile(exp(ind_llh), 0.50); % Median likelihood value;
    eps_range = 0:eps_max/step:eps_max;
    
    intl_lambda = hlambda_1 + randn(p,k);
    intl_psi = hpsi_1;
  
    for iter = 1:maxiter
    
        % Get REM estimates;
        [hlambda_2, hpsi_2, est_gamma(iter), weights] = RobustEMAlg(X,k,eps_range(iter),intl_lambda,intl_psi);
                        
        % Calculate prob. that weights <= 1/2;
        chk(iter) = checkEpsFA(hlambda_2,hpsi_2,est_gamma(iter),eps_range(iter));
        
        if chk(iter) > delta 
            break;
        end
    end
    
    % Truncate gamma and chk;
    chk = chk(1:iter);
    est_gamma = est_gamma(1:iter);
           
    % Compute and store congruence coefficient
    R_classic(l) = computeCongruence(hlambda_1, lambda_t); % first matrix gets procrustes rotation;
    R_robust(l) = computeCongruence(hlambda_2, lambda_t);
    
    R_classic_b(l) = computeCongruence(hlambda_1, lambda_b); % first matrix gets procrustes rotation;
    R_robust_b(l) = computeCongruence(hlambda_2, lambda_b);
            
    % Store gamma value
    gamma_values(l) = est_gamma(iter);  
    
    % Store epsilon choice
    opt_eps(l) = eps_range(iter);

end


%%%%%%%%%%% Data Aggregation %%%%%%%%%%%%%%%%%%%%%%

meanR = mean(R_classic);
stdR = std(R_classic,1,1);
meanR_robust = mean(R_robust);
stdR_robust = std(R_robust,1,1);

meanR_b = mean(R_classic_b);
stdR_b = std(R_classic_b,1,1);
meanR_robust_b = mean(R_robust_b);
stdR_robust_b = std(R_robust_b,1,1);

meanGamma = mean(gamma_values);
stdGamma = std(gamma_values,1,1);

R_EM = [meanR stdR];
R_REM = [meanR_robust stdR_robust];
gamma = [meanGamma stdGamma];

R_EM_minor = [meanR_b stdR_b];
R_REM_minor = [meanR_robust_b stdR_robust_b];
   
end