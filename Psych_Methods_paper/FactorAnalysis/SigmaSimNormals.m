function [sigma, lambda_common, lambda_unique] = SigmaSimNormals(p,k,communality)
% p = # variables (items)
% k = # factors

Z = randn(p,k);

g = 1./vecnorm(Z,2,2);

% Scale to set communality
if communality == 1
    B1 = diag(1/10*randi([2,4],p,1));
elseif communality == 2
    B1 = diag(1/10*randi([2,8],p,1));
elseif communality == 3
    B1 = diag(1/10*randi([6,8],p,1));
else 
    B1 = zeros(p,p);
end

B2 = eye(p,p) - B1;

% Final factor loading matrix
lambda_common = sqrt(B1)*g.*Z;
lambda_unique = sqrt(B2);

sigma = lambda_common*lambda_common' + lambda_unique*lambda_unique';

end