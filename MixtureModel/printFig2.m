function printFig2()

clear; close all;

set(0,'defaulttextinterpreter','latex')

% Example 1;
load('GMM_output_Sim_2_k_2_delta_1.mat','weights'); 
 
subplot(2,2,1)
histogram(weights(weights(:,2) == 0,1),10, ...
    'Normalization','probability','FaceColor','#A2142F')
xlim([0,1]); xticks(0:0.25:1); ylim([0,1]); yticks(0:0.25:1);
title([{'$$\textbf{Example 1}$$'},{''}]);
ylabel({'$$\textbf{Majority} $$';'';'';'Frequency'}); xlabel([]);
PrettyFig
legend off

subplot(2,2,3)
histogram(weights(weights(:,2) == 1,1),10, ...
    'Normalization','probability','FaceColor','#EDB120')
ylim([0,1]); yticks(0:0.25:1); xlim([0,1]); xticks(0:0.25:1); 
ylabel(({'$$\textbf{Minority} $$';'';'';'Frequency'})); xlabel('Weight');
PrettyFig
legend off


% Example 2;
load('GMM_output_Sim_3_k_2_delta_1.mat','weights'); 

subplot(2,2,2)
histogram(weights(weights(:,2) == 0,1),10, ...
    'Normalization','probability','FaceColor','#A2142F')
ylim([0,1]); yticks(0:0.25:1); xlim([0,1]); xticks(0:0.25:1); 
ylabel([]); xlabel([]);
title([{'$$\textbf{Example 2}$$'},{''}]);
PrettyFig
legend off

subplot(2,2,4)
histogram(weights(weights(:,2) == 1,1),10, ...
    'Normalization','probability','FaceColor','#EDB120')
ylim([0,1]); yticks(0:0.25:1); xlim([0,1]); xticks(0:0.25:1); 
ylabel([]); xlabel('Weight');
PrettyFig
legend off


set(gcf,'Position',[50 50 900 800])

print('Fig2','-depsc','-r600','-painters')
%print('Fig2','-dpng','-r600','-painters')

end