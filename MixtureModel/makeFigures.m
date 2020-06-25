function [est_gmm] = makeFigures(X,est_gmm)

k = est_gmm.NumComponents;
[~, rnk] = sort(est_gmm.ComponentProportion);
    
th = 0:2*pi/50:2*pi;
x = cos(th); y = sin(th);

hold on
for i = 1:k
    ell = 2*chol(est_gmm.Sigma(:,:,rnk(i)))'*[x; y] + est_gmm.mu(rnk(i),:)';
    fill(ell(1,:), ell(2,:), [.7,0.1,1],'FaceAlpha',0.2,'LineWidth',1.5)
end
scatter(X(1,:), X(2,:),25,[0.45,0.45,0.45],'Marker','.')
scatter(est_gmm.mu(:,1), est_gmm.mu(:,2),40,[0,0,0],'Marker','x','LineWidth',2.5)

PrettyFig
legend off
hold off

end