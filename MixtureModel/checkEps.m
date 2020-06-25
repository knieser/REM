function chk = checkEps(mu,sigma,mix,gamma,eps)
%{
This function checks whether the selected epsilon value leads to REM estimates 
that meet the heuristic for hyperparameter selection.
    
INPUT:
    mu: (p x k) estimated mean vectors
    sigma: (p x p x k) estimated covariance matrices
    mix: (1 x p) estimated mixture proportions
    gamma: estimated gamma from REM estimation
    eps: hyperparameter for REM estimation 

OUTPUT:
    chk: probability that gamma.*f <= (1-gamma)*eps 
         or equivalently weights <= 1/2
%}  
     
% Generate data;
n = 10000;
p = size(mu,1);
k = size(mu,2);
gm = gmdistribution(mu',sigma,mix);
X = random(gm,n);
X = X';
    
% Calculate likelihood
for j = 1:k
    log_g(:,j) = -1/2*(p*log(2*pi) + log(det(sigma(:,:,j))) + diag((X - mu(:,j))'/sigma(:,:,j)*(X - mu(:,j))));   
end
f = sum(mix.*exp(log_g),2); 

% Calculate weight distribution;
weights = gamma*f ./ (gamma.*f + (1-gamma)*eps);

% Calculate chk;
chk = mean(weights < 1/2);
    
end