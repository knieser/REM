function [RMSE_mu, RMSE_sigma] = MSE(gmm, gmm_ref)
%{
This function takes two GMM distributions and calculates the RMSE between
the means and the Frobenius norm of the difference between the covariance
matrices.
%}

[~,idx0] = sort(gmm_ref.ComponentProportion);
[~,idx1] = sort(gmm.ComponentProportion);

diff_mu = gmm.mu(idx1,:) - gmm_ref.mu(idx0,:);
RMSE_mu = sqrt(mean((diff_mu(:)).^2));

diff_sigma = gmm.Sigma(:,:,idx1) - gmm_ref.Sigma(:,:,idx0);
RMSE_sigma = sqrt(mean((reshape(diff_sigma,1,[]).^2)));

end