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
mix = [0.6 0.3 0.1];
   

%%%%%%% Make data %%%%%%%%%

[X, grp_flag] = LPA_data_sim(sim,p_0,k_0,n,mix);


%%%%%%% Run simulation %%%%%%%

[est_gmm, est_gmm_2, est_weights, est_gamma, opt_eps] = LPA_estimates(X,k,delta,global_iter,step);

% Add flag for minority group;
minority = [grp_flag grp_flag > k];
weights = [est_weights minority(:,2)];


%%%%% Make figures %%%%%%

figure;
hold on
scatter(X(1,:), X(2,:),50,[0.3,0.3,0.3],'Marker','.')
xlabel('Domain A')
ylabel('Domain B')
PrettyFig
legend off
title('Simulated Data')
xlim([0,10])
ylim([0,10])
hold off
%print('Sim3_Data','-dpng','-r300');


makeFigures(X,est_gmm);
title('Classic')
xlim([0,10])
ylim([0,10])
%print('Sim1_Classic_K3','-dpng','-r300')


makeFigures(X,est_gmm_2);
title('Robust')
xlim([0,10])
ylim([0,10])
%print('Sim1_Robust_K3','-dpng','-r300')

figure;
hold on
histogram(weights(weights(:,2) == 0,1))
histogram(weights(weights(:,2) == 1,1))
xlabel('Weight')
ylabel('Count')
PrettyFig
legend('Majority','Minority','Location','northwest')
title('Weight Distribution')

    
save LPA_output_5_15_20;

end