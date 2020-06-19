function chk = checkEpsFA(lambda,psi,gamma,eps)

n = 10000;

p = size(lambda,1);

% Generate some data;
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