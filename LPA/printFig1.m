function printFig1()

clear; close all;

load('LPA_output_Sim_2_k_2_delta_1.mat'); 
 
subplot(5,2,1)
makeFigures(X,true_gmm);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);
title('Example 1','FontWeight','bold')

subplot(5,2,3)
makeFigures(X,est_gmm);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);

subplot(5,2,5)
makeFigures(X,est_gmm_2);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);

subplot(5,2,7)
histogram(weights(weights(:,2) == 0,1),10, ...
    'Normalization','probability','FaceColor','#A2142F')
xlabel('Weight'); xlim([0,1]); xticks(0:0.25:1); 
ylabel('Frequency'); ylim([0,1]); yticks(0:0.25:1);
PrettyFig
legend off

subplot(5,2,9)
histogram(weights(weights(:,2) == 1,1),10, ...
    'Normalization','probability','FaceColor','#EDB120')
xlabel('Weight'); xlim([0,1]); xticks(0:0.25:1); 
ylabel('Frequency'); ylim([0,1]); yticks(0:0.25:1);
PrettyFig
legend off

load('LPA_output_Sim_3_k_2_delta_1.mat'); 
 
subplot(5,2,2)
makeFigures(X,true_gmm);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);
title('Example 2','FontWeight','bold')

subplot(5,2,4)
makeFigures(X,est_gmm);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);

subplot(5,2,6)
makeFigures(X,est_gmm_2);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);

subplot(5,2,8)
histogram(weights(weights(:,2) == 0,1),10, ...
    'Normalization','probability','FaceColor','#A2142F')
xlabel('Weight'); xlim([0,1]); xticks(0:0.25:1); 
ylabel('Frequency'); ylim([0,1]); yticks(0:0.25:1);
PrettyFig
legend off

subplot(5,2,10)
histogram(weights(weights(:,2) == 1,1),10, ...
    'Normalization','probability','FaceColor','#EDB120')
xlabel('Weight'); xlim([0,1]); xticks(0:0.25:1); 
ylabel('Frequency'); ylim([0,1]); yticks(0:0.25:1);
PrettyFig
legend off

set(gcf,'Position',[50 50 600 1000])

print('Fig1','-depsc','-r600','-painters')

end