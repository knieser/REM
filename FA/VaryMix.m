function [R_EM, R_REM, gamma, opt_eps, R_EM_minor, R_REM_minor] = VaryMix(p_0,k_0,n,k_choice,corrupt_pct,sparse,communality,m,delta,step,n_sim,seed)
    
% Fix k;
k = k_0;

% Initialize;
R_EM = zeros(length(corrupt_pct), 2);
R_REM = zeros(length(corrupt_pct), 2);

R_EM_minor = zeros(length(corrupt_pct), 2);
R_REM_minor = zeros(length(corrupt_pct), 2);

gamma = zeros(length(corrupt_pct), 2);
opt_eps = zeros(length(corrupt_pct), n_sim);

% Vary corruption;
for c = 1:length(corrupt_pct)
    
    msg = ['Working on corrupt pct = ', num2str(100*corrupt_pct(c)),'%'];
    disp(msg)
    
    [R_EM(c,:), R_REM(c,:), gamma(c,:), opt_eps(c,:), R_EM_minor(c,:), R_REM_minor(c,:)] = ....
    FA_Hetero_Pop_Sim(p_0,k_0,n,k,corrupt_pct(c),sparse,communality,m,delta,step,n_sim,seed);
end


%%%% Outputs %%%%

% R_(R)EM(:,1): mean R over n_sim simulations;
% R_(R)EM(:,2): SE of R over n_sim simulations;
% gamma = [mean SE];


%%%%%%%%%%% Make Figures %%%%%%%%%%%%%%%%%%%

maj_pct = 1 - corrupt_pct;

% Plot congruence vs mix %;
figure;
hold on
errorbar(maj_pct, R_EM(:,1), R_EM(:,2),'LineWidth',1.5,'Color','b')
errorbar(maj_pct, R_REM(:,1), R_REM(:,2),'LineWidth',1.5,'Color','r')
set(gca,'YLim',[0 1])  
ylabel('Congruence Coeff')
xlabel('Majority %')
PrettyFig
legend('Classic','Robust', 'Location','southeast','Orientation','vertical')
hold off
%print('Heterogeneity high communality m0','-dpng','-r300')


% Plot congruence to minority vs mix %;
figure;
hold on
errorbar(maj_pct, R_EM_minor(:,1), R_EM_minor(:,2),'LineWidth',1.5,'Color','b')
errorbar(maj_pct, R_REM_minor(:,1), R_REM_minor(:,2),'LineWidth',1.5,'Color','r')
set(gca,'YLim',[0 1])  
ylabel('Congruence Coeff')
xlabel('Majority %')
PrettyFig
legend('Classic','Robust', 'Location','southeast','Orientation','vertical')
hold off
%print('Heterogeneity high communality m0','-dpng','-r300')


% Plot estimated gamma vs mix %;
figure
hold on
errorbar(maj_pct, gamma(:,1), gamma(:,2),'LineWidth',1.5,'Color','m')
plot(0:0.01:1, 0:0.01:1,'--','LineWidth',1,'Color','k')
hold off
ylabel('Estimated Gamma')
xlabel('Majority %')
set(gca,'YLim',[0.6 1])  
PrettyFig
legend('off')
%print('Gamma ep16 n2000','-dpng','-r300')

  
end