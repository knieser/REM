function printFig2()

clear; close all;

load('LPA_output_Sim_1_k_1_delta_1.mat'); 
 
subplot(3,3,1)
makeFigures(X,true_gmm);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);
title('K = 1','FontWeight','bold')

subplot(3,3,4)
makeFigures(X,est_gmm);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);

subplot(3,3,7)
makeFigures(X,est_gmm_2);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);

load('LPA_output_Sim_1_k_2_delta_1.mat'); 
 
subplot(3,3,2)
makeFigures(X,true_gmm);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);
title('K = 2','FontWeight','bold')

subplot(3,3,5)
makeFigures(X,est_gmm);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);

subplot(3,3,8)
makeFigures(X,est_gmm_2);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);


load('LPA_output_Sim_1_k_3_delta_1.mat'); 
 
subplot(3,3,3)
makeFigures(X,true_gmm);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);
title('K = 3','FontWeight','bold')

subplot(3,3,6)
makeFigures(X,est_gmm);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);

subplot(3,3,9)
makeFigures(X,est_gmm_2);
xlim([0,10]); xticks(0:5:10); ylim([0,10]); yticks(0:5:10);

set(gcf,'Position',[50 50 1000 1000])

print('Fig2','-depsc','-r600','-painters')

end