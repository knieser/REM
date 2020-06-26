function chk = checkEpsFA(lambda,psi,gamma,eps)
%{
This function checks whether the selected epsilon value leads to REM estimates 
that meet the heuristic for hyperparameter selection.
    
INPUT:
    lambda: (p x k) estimated lambda
    psi: (p x 1) estimated psi
    gamma: estimated gamma from REM estimation
    eps: hyperparameter for REM estimation 

OUTPUT:
    chk: probability that gamma.*f <= (1-gamma)*eps 
         or equivalently weights <= 1/2
%}


% Generate some data;
n = 10000;
p = size(lambda,1);
sigma = lambda*lambda' + diag(psi);
V = chol(sigma);
X = randn(n,p)*V;
X = X';

% Calculate likelihood;
logdetV = log(det(V));  % the det(V) is the square root of the det of sigmax;
log_f = -(1/2)*(p*log(2*pi)+2*logdetV+diag((X'/V)/(V')*X));
f = exp(log_f);

% Calculate weight distribution;
weights = gamma*f ./ (gamma*f + (1-gamma)*eps);

% Calculate chk;
chk = mean(weights < 1/2);

    
end