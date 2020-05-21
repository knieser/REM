function printFig4()

clear; close all;

load('FA_output_Sim_2_k_5_comm_3_delta_5.mat','k_choice','R_EM','R_REM'); 

hold on
errorbar(k_choice, R_EM(:,1), R_EM(:,2),'LineWidth',1.5,'Color','b')
errorbar(k_choice, R_REM(:,1), R_REM(:,2),'LineWidth',1.5,'Color','r')
set(gca,'YLim',[0.5 1])  
ylabel('RV Coefficient')
xlabel('K Choice')
PrettyFig
legend('EM','REM', 'Location','southeast','Orientation','horizontal')
hold off

%set(gcf,'Position',[50 50 600 800])

print('Fig4','-depsc','-r600','-painters')

end