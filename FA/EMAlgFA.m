function [hlambda, hpsi, ind_llh] = EMAlg(X,k)

p = size(X,1);
n = size(X,2);

Cxx = (1/n)*(X*X'); 

% Initial guess of lambda and psi;
[hlambda, hpsi] = factoran(X',k, 'maxit', 1000);

% Tolerance parameters for EM alg
tol = 1e-15;
maxiter = 500;
llh = -inf(1,maxiter);

for iter = 2:maxiter
    
    % log-likelihood
    V = chol(diag(hpsi)+hlambda*hlambda');  % diag(hpsi) + hlambda*hlambda' = V'V;
    logdetV = log(det(V));                  % the det(V) is the square of the det of sigmax;
    loglike = -(1/2)*(p*log(2*pi)+2*logdetV+diag((X'/V)/(V')*X));
     
    llh(iter) = -(n/2)*(p*log(2*pi)+2*logdetV+(1/n)*trace((X'/V)/(V')*X));
    if abs(llh(iter) - llh(iter-1)) < tol*abs(llh(iter-1)); break; end
        
    % EM algorithm
    hbeta   = (hlambda'/V)/V';
    chol_Ezz = chol(eye(k) - hbeta*hlambda + hbeta*Cxx*hbeta');
    
    hlambda = ((Cxx*hbeta')/chol_Ezz)/chol_Ezz';
    hpsi    = diag((eye(p) - hlambda*hbeta)*Cxx);
    
end

ind_llh = loglike; %log likelihood for each individual

end