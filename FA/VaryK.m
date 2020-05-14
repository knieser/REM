function [R_EM, R_REM, gamma, opt_eps, R_EM_minor, R_REM_minor] = VaryK(p_0,k_0,n,k_choice,corrupt_pct,sparse,communality,m,delta,step,n_sim,seed)
   
% Fix corruption;
corrupt_pct = 0.10;

% Initialize;
R_EM = zeros(length(k_choice), 2);
R_REM = zeros(length(k_choice), 2);

R_EM_minor = zeros(length(k_choice), 2);
R_REM_minor = zeros(length(k_choice), 2);

gamma = zeros(length(k_choice), 2);
opt_eps = zeros(length(k_choice), n_sim);

% Vary choice of K;
for k = 1:length(k_choice)
    
    msg = ['Working on k = ', num2str(k_choice(k))];
    disp(msg)
    
    [R_EM(k,:), R_REM(k,:), gamma(k,:), opt_eps(k,:), R_EM_minor(k,:), R_REM_minor(k,:)] = ....
    FA_Hetero_Pop_Sim(p_0,k_0,n,k_choice(k),corrupt_pct,sparse,communality,m,delta,step,n_sim,seed);
end

%%%% Outputs %%%%

% R_(R)EM(:,1): mean R over n_sim simulations;
% R_(R)EM(:,2): SE of R over n_sim simulations;
% gamma = [mean SE];


%%%%%%%%%%% Make Figures %%%%%%%%%%%%%%%%%%%

% Plot congruence vs mix %;
figure;
hold on
errorbar(k_choice, R_EM(:,1), R_EM(:,2),'LineWidth',1.5,'Color','b')
errorbar(k_choice, R_REM(:,1), R_REM(:,2),'LineWidth',1.5,'Color','r')
set(gca,'YLim',[0 1])  
ylabel('Congruence Coeff')
xlabel('K Choice')
PrettyFig
legend('Classic','Robust', 'Location','southeast','Orientation','vertical')
hold off
%print('Heterogeneity high communality m0','-dpng','-r300')


% Plot estimated gamma vs mix %;
figure
hold on
errorbar(k_choice, gamma(:,1), gamma(:,2),'LineWidth',1.5,'Color','m')
hold off
ylabel('Estimated Gamma')
xlabel('K Choice')
set(gca,'YLim',[0 1])  
PrettyFig
legend('off')
%print('Gamma ep16 n2000','-dpng','-r300')



    
end