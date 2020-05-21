function [est_gmm] = makeFigures(X,est_gmm)

th = 0:2*pi/50:2*pi;
x = cos(th); y = sin(th);

hold on
scatter(X(1,:), X(2,:),25,[0.45,0.45,0.45],'Marker','.')

for i = 1:size(est_gmm.Sigma,3)
    ell = 2*chol(est_gmm.Sigma(:,:,i))'*[x; y] + est_gmm.mu(i,:)';
    fill(ell(1,:), ell(2,:), [.7,0.1,1],'FaceAlpha',0.2,'LineWidth',1.5)
end
scatter(est_gmm.mu(:,1), est_gmm.mu(:,2),40,[0,0,0],'Marker','x','LineWidth',2.5)

xlabel('Domain A')
ylabel('Domain B')
PrettyFig
legend off
hold off

end