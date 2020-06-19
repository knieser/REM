function printFig3()

clear; close all;

set(0,'defaulttextinterpreter','latex')

% K = 1
load('LPA_output_Sim_1_k_1_delta_1.mat','X','true_gmm','est_gmm','est_gmm_2'); 
 
subplot(3,3,1)
makeFigures(X,true_gmm);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);
title([{'$$\textbf{K = 1}$$'},{''}]);
ylabel({'$$\textbf{Simulated}$$';'';'';'Domain A'}); xlabel([]);

subplot(3,3,4)
makeFigures(X,est_gmm);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);
ylabel({'$$\textbf{EM}$$';'';'';'Domain A'}); xlabel([]);

subplot(3,3,7)
makeFigures(X,est_gmm_2);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);
ylabel({'$$\textbf{REM} $$';'';'';'Domain A'}); xlabel('Domain B');


% K = 2
load('LPA_output_Sim_1_k_2_delta_1.mat','X','true_gmm','est_gmm','est_gmm_2');  
 
subplot(3,3,2)
makeFigures(X,true_gmm);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);
title([{'$$\textbf{K = 2}$$'},{''}]);
ylabel([]); xlabel([]);

subplot(3,3,5)
makeFigures(X,est_gmm);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);
ylabel([]); xlabel([]);

subplot(3,3,8)
makeFigures(X,est_gmm_2);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);
ylabel([]); xlabel('Domain B');


% K = 3
load('LPA_output_Sim_1_k_3_delta_1.mat','X','true_gmm','est_gmm','est_gmm_2');
 
subplot(3,3,3)
makeFigures(X,true_gmm);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);
title([{'$$\textbf{K = 3}$$'},{''}]);
ylabel([]); xlabel([]);

subplot(3,3,6)
makeFigures(X,est_gmm);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);
ylabel([]); xlabel([]);

subplot(3,3,9)
makeFigures(X,est_gmm_2);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);
ylabel([]); xlabel('Domain B');


set(gcf,'Position',[50 50 1000 1000])

print('Fig3','-depsc','-r600','-painters')
%print('Fig3','-dpng','-r600','-painters')

end