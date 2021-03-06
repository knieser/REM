function printFig1(delta)

set(0,'defaulttextinterpreter','latex')

% Example 1;
load(['GMM_output_Sim_2_k_2_delta_',num2str(100*delta),'.mat'],'X','true_gmm','est_gmm','est_gmm_2'); 
 
subplot(3,2,1)
makeFigures(X,true_gmm);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);
title([{'$$\textbf{Example 1}$$'},{''}]);
ylabel({'$$\textbf{Simulated} $$';'';'';'Domain B'}); xlabel([]);

subplot(3,2,3)
makeFigures(X,est_gmm);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);
ylabel({'$$\textbf{EM}$$';'';'';'Domain B'}); xlabel([]);

subplot(3,2,5)
makeFigures(X,est_gmm_2);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);
ylabel({'$$\textbf{REM}$$';'';'';'Domain B'}); xlabel('Domain A');


% Example 2
load(['GMM_output_Sim_3_k_2_delta_',num2str(100*delta),'.mat'],'X','true_gmm','est_gmm','est_gmm_2'); 
 
subplot(3,2,2)
makeFigures(X,true_gmm);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);
title([{'$$\textbf{Example 2}$$'},{''}]);
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

print('Fig1','-depsc','-r600','-painters')

end