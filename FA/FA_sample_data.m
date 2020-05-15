function [X,lambda_target_1,lambda_target_2] = FA_sample_data(sigma_01, sigma_02, lambda_01, lambda_02,n,corrupt_pct)

p = size(sigma_01,1);
corrupt_number = ceil(n*(corrupt_pct));

% Sample data from multivariate normal with mean 0 and covariance sigma;
R_01 = chol(sigma_01);
X_01 = randn(n,p)*R_01;
X_01 = X_01';

R_02 = chol(sigma_02);
X_02 = randn(n,p)*R_02;
X_02 = X_02';

X = [X_01(:,1:n-corrupt_number) X_02(:,1:corrupt_number)];  
V = var(X,0,2);

X = normalize(X')';

% Scale lambda_01 to align with normalized data;
invV = 1./V;
scaled_lambda_01 = sqrt(invV).*lambda_01;
scaled_lambda_02 = sqrt(invV).*lambda_02;

 % Just focus on the "main" lambda at this point (ie lambda_01);
lambda_target_1 = scaled_lambda_01;  
lambda_target_2 = scaled_lambda_02;
    
end