function [R] = computeCongruence(lambda_1, lambda_2)

% Calculate "squares" of matrices XX';
sq_lambda_1 = lambda_1*lambda_1';
sq_lambda_2 = lambda_2*lambda_2';

% Calculate RV coeffs;
R = trace(sq_lambda_1*sq_lambda_2)/sqrt(...
    trace(sq_lambda_1*sq_lambda_1')*trace(sq_lambda_2*sq_lambda_2')...
    );
        
end