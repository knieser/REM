function [R] = computeCongruence(lambda_1, lambda_2)

% RV Coeff
% Calculate "squares" of matrices XX';
sq_lambda_1 = lambda_1*lambda_1';
sq_lambda_2 = lambda_2*lambda_2';

% Calculate RV coeffs;
R = trace(sq_lambda_1*sq_lambda_2)/sqrt(...
    trace(sq_lambda_1*sq_lambda_1')*trace(sq_lambda_2*sq_lambda_2')...
    );
    
%{   
% Make sure lambda_1 and lambda_2 are aligned;
if size(lambda_1,2) == 1
    if lambda_1(1)*lambda_2(1) < 0
        lambda_1 = -1*lambda_1; 
    end
else
    lambda_1 = rotatefactors(lambda_1,'method','procrustes','target',lambda_2,'type','oblique');
end
    
% Tucker Congruence
R = trace(lambda_1*lambda_2') / sqrt(trace(lambda_1*lambda_1')*trace(lambda_2*lambda_2'));  
  
%}
        
end