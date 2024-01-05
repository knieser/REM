function [hlambda, hpsi, gamma, weights, alt_llh, llh] = RobustEMAlg(X,k,epsilon,intl_lambda,intl_psi)
%{
This function obtains REM estimates.
    
INPUT:
    X: (p x n) dataset
    k: number of latent factors
    epsilon: hyperparameter
    intl_lambda: (p x k) intial value for lambda
    intl_psi: (p x 1) initial value for psi
    
OUTPUT:
    hlambda: (p x k) estimated lambda
    psi: (p x 1) estimated psi
    gamma: REM estimated gamma
    weights: individual-level probabilistic weights
    alt_llh: (n x 1) individual-level modified loglikelihood values 
%} 

% Get data dimensions;
p = size(X,1);
n = size(X,2);

% Initial guesses;
hlambda = intl_lambda;
hpsi = intl_psi;
gamma = 0.5;

% Tolerance parameters
tol = 1e-15;
maxiter = 500;
llh = -inf(1,maxiter);

for iter = 2:maxiter
    
    try 
        chol(diag(hpsi)+hlambda*hlambda');      
    catch 
        warning(['Set gamma to 0 for log(epsilon) = ',num2str(log(epsilon))])
        gamma = 0;
        weights = zeros(n,1);
        break;
    end
    
    V = chol(diag(hpsi)+hlambda*hlambda'); % diag(hpsi) + hlambda*hlambda' = V'V;
    logdetV = log(det(V));  % the det(V) is the square root of the det of sigmax;
    
    % Log-likelihood
    ind_llh = -(1/2)*(p*log(2*pi)+2*logdetV+diag((X'/V)/(V')*X));
    alt_ind_llh = log(gamma*exp(ind_llh) + (1-gamma)*epsilon);
    llh(iter) = sum(alt_ind_llh);   
    
    if abs(llh(iter) - llh(iter-1)) < tol*abs(llh(iter-1)); break; end
    
    % Create weights
    scaledllh = log(gamma)+ind_llh; % Add in gamma

    weights   = log(1-gamma) + log(epsilon) - scaledllh;
    weights   = 1./( 1 + exp(weights) );

    gamma = 1/n*sum(weights);
    Cxx = (X*(weights.*X'))/(n*gamma);
    
    hbeta   = (hlambda'/V)/V';
    chol_Ezz = chol(eye(k) - hbeta*hlambda + hbeta*Cxx*hbeta');
    
    hlambda = ((Cxx*hbeta')/chol_Ezz)/chol_Ezz';
    hpsi    = diag((eye(p) - hlambda*hbeta)*Cxx);

end

% Store alternate and traditional log-likelihood value
alt_llh = sum(alt_ind_llh);
llh     = sum(ind_llh);

end