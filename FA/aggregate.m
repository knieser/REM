function [A_agg] = aggregate(A)

meanA = mean(A);
stdA = std(A,1,1);

A_agg = [meanA stdA]';
  
end