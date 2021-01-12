function printFig7(delta)

set(0,'defaulttextinterpreter','latex')

% Example 5;
load(['GMM_output_Sim_6_k_2_delta_',num2str(100*delta),'.mat'],'X','true_gmm','est_gmm','est_gmm_2'); 
 
subplot(3,1,1)
makeFigures(X,true_gmm);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);
title([{'$$\textbf{Example 5}$$'},{''}]);
ylabel({'$$\textbf{Simulated} $$';'';'';'Domain B'}); xlabel([]);

subplot(3,1,2)
makeFigures(X,est_gmm);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);
ylabel({'$$\textbf{EM}$$';'';'';'Domain B'}); xlabel([]);

subplot(3,1,3)
makeFigures(X,est_gmm_2);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);
ylabel({'$$\textbf{REM}$$';'';'';'Domain B'}); xlabel('Domain A');

set(gcf,'Position',[50 50 500 800])

print('Fig7','-depsc','-r600','-painters')

end