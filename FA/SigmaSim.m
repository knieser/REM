function [sigma, lambda_common, lambda_unique] = SigmaSim(p,k,communality)
% p = # variables (items)
% k = # factors

% Need row entries to sum to k-1;
A = zeros(p,k);

% First choose randomly the values for first column
A(:,1) = randi([0,k-1],p,1);

% Go through each entry by row
for i = 1:p 
    for j = 2:k-1
        if sum(A(i,:),2) == k-1 
            A(i,j) = 0;
        else
            A(i,j) = randi([0,(k-1) - sum(A(i,:),2)]);
        end
    end
end      

% Set last column to be remaining balance
A(:,k) = (k-1) - sum(A,2);

% Add normal deviation
c = 1/10*randi([7,9],k,1); % c = 0.7, 0.8, 0.9 with equal prob
x1 = randn(p,k);
d = 1./vecnorm(x1,2,2);

Y = A.*c' + d.*x1.*sqrt(1-c.^2)';

% Apply skewing function
Y2 = Y + abs(Y) + 0.2;
Y3 = abs(Y)+ 0.2;
Z = (1.2/2.2) * (Y.*Y2) ./ Y3;

g = 1./vecnorm(Z,2,2);

% Scale to set communality
if communality == 1
    B1 = diag(1/10*randi([2,4],p,1));
elseif communality == 2
    B1 = diag(1/10*randi([2,8],p,1));
elseif communality == 3
    B1 = diag(1/10*randi([6,8],p,1));
elseif communality == 4
    B1 = diag(1/10*randi([9,10],p,1));
else 
    B1 = zeros(p,p);
end

B2 = eye(p,p) - B1;

% Final factor loading matrix
lambda_common = sqrt(B1)*g.*Z;
lambda_unique = sqrt(B2);

sigma = lambda_common*lambda_common' + lambda_unique*lambda_unique';

end