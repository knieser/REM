function [R] = computeCongruence(lambda_1, lambda_2)

% Calculate "squares" of matrices XX';
sq_lambda_1 = lambda_1*lambda_1';
sq_lambda_2 = lambda_2*lambda_2';

% Calculate RV coeffs;
R = trace(sq_lambda_1*sq_lambda_2)/sqrt(...
    trace(sq_lambda_1*sq_lambda_1')*trace(sq_lambda_2*sq_lambda_2')...
    );

    
%{
nvar = size(lambda_1,1);
nfactor = size(lambda_1,2);

rot_lambda_1 = rotatefactors(lambda_1,'Method','procrustes','target',lambda_2); 
rot_lambda_2 = lambda_2;
%rot_lambda_1 = rotatefactors(lambda_1,'Method','quartimax');
%rot_lambda_2 = rotatefactors(lambda_2,'Method','quartimax');


% Tucker congruence coefficient
Rk = 1:nfactor;

for k = 1:length(Rk)
    Rk(k) = (rot_lambda_1(:,k)'*rot_lambda_2(:,k)) / ...
        sqrt(rot_lambda_1(:,k)'*rot_lambda_1(:,k)*rot_lambda_2(:,k)'*rot_lambda_2(:,k));
end

%R = mean(Rk);

%RMSE
bias = 0;
for p = 1:nvar
    for k = 1:nfactor
        bias = bias + rot_lambda_1(p,k) - rot_lambda_2(p,k);
    end
end
%}

         
end