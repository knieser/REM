function [est_gmm,est_gmm_2, weights, est_gamma, opt_eps] = LPA_sim_main(sim,k,delta)
% Sim Key
% 1 = No corruption
% 2 = Scattered minority
% 3 = Scattered minority w cross

% e.g. LPA_sim_main(1,2,0.05);
  
%%%%%%%% Parameters %%%%%%%%%%%

% Number of starting points for global optimization;
global_iter = 300;

% Number of iterations to search for epsilon;
step = 100;

% Simulated data parameters;
p_0 = 2; % number of dimensions of data;
k_0 = 3; % true number of clusters;
n = 1000;
mix = [0.70 0.20 0.10];
   

%%%%%%% Make data %%%%%%%%%

[X, true_gmm, grp_flag] = LPA_data_sim(sim,p_0,k_0,n,mix);


%%%%%%% Run simulation %%%%%%%

[est_gmm, est_gmm_2, est_weights, est_gamma, opt_eps] = LPA_estimates(X,k,delta,global_iter,step);

% Add flag for minority group;
minority = [grp_flag grp_flag > k];
weights = [est_weights minority(:,2)];


%%%%% Make figures %%%%%%

true_gmm = gmdistribution(true_gmm.mu(1:k,:), true_gmm.Sigma(:,:,1:k), true_gmm.ComponentProportion(1:k));

makeFigures(X,true_gmm);
%title('Simulated Data')
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);
print1 = ['Sim_',num2str(sim),'_k_',num2str(k),'_delta_', num2str(100*delta),'pct_true_gmm_',date];
print(print1,'-dpng','-r300');

makeFigures(X,est_gmm);
%title('Classic')
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);
print2 = ['Sim_', num2str(sim),'_k_',num2str(k),'_delta_', num2str(100*delta),'pct_EM_gmm_',date];
print(print2,'-dpng','-r300');

makeFigures(X,est_gmm_2);
%title('Robust')
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);
print3 = ['Sim_', num2str(sim),'_k_',num2str(k),'_delta_', num2str(100*delta),'pct_REM_gmm_',date];
print(print3,'-dpng','-r300');

figure;
hold on
histogram(weights(weights(:,2) == 0,1),20,'Normalization','probability')
histogram(weights(weights(:,2) == 1,1),20,'Normalization','probability')
xlabel('Weight'); xlim([0,1]); xticks(0:0.25:1); 
ylabel('Frequency'); ylim([0,1]); yticks(0:0.25:1);
PrettyFig
legend('Majority','Minority','Location','southoutside','Orientation','horizontal','FontSize',10)
%title('Weight Distribution')
print4 = ['Sim_', num2str(sim),'_k_',num2str(k),'_delta_', num2str(100*delta),'pct_weights_',date];
print(print4,'-dpng','-r300');

close all;

print5 = ['LPA_output_Sim_', num2str(sim),'_k_',num2str(k),'_delta_', num2str(100*delta)'];
save(print5);

end