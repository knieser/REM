function printFig6(delta)

set(0,'defaulttextinterpreter','latex')

% Example 3;
load(['GMM_output_Sim_4_k_2_delta_',num2str(100*delta),'.mat'],'X','true_gmm','est_gmm','est_gmm_2'); 
 
subplot(3,2,1)
makeFigures(X,true_gmm);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);
title([{'$$\textbf{Example 3}$$'},{''}]);
ylabel({'$$\textbf{Simulated} $$';'';'';'Domain B'}); xlabel([]);

subplot(3,2,3)
makeFigures(X,est_gmm);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);
ylabel({'$$\textbf{EM}$$';'';'';'Domain B'}); xlabel([]);

subplot(3,2,5)
makeFigures(X,est_gmm_2);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);
ylabel({'$$\textbf{REM}$$';'';'';'Domain B'}); xlabel('Domain A');


% Example 4
load(['GMM_output_Sim_5_k_2_delta_',num2str(100*delta),'.mat'],'X','true_gmm','est_gmm','est_gmm_2'); 
 
subplot(3,2,2)
makeFigures(X,true_gmm);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);
title([{'$$\textbf{Example 4}$$'},{''}]);
ylabel([]); xlabel([]);

subplot(3,2,4)
makeFigures(X,est_gmm);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);
ylabel([]); xlabel([]);

subplot(3,2,6)
makeFigures(X,est_gmm_2);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);
ylabel([]); xlabel('Domain A');

set(gcf,'Position',[50 50 800 1000])

print('Fig6','-depsc','-r600','-painters')

end