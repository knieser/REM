function printFig3()

clear; close all;

load('FA_output_Sim_1_k_1_comm_3_delta_5.mat','corrupt_pct','R_EM','R_REM','gamma_values'); 
 
subplot(3,2,1)
makeCongruenceFig(corrupt_pct,R_EM,R_REM)
     
subplot(3,2,2)
makeGammaFig(corrupt_pct,gamma_values) 

load('FA_output_Sim_1_k_1_comm_2_delta_5.mat','corrupt_pct','R_EM','R_REM','gamma_values'); 
 
subplot(3,2,3)
makeCongruenceFig(corrupt_pct,R_EM,R_REM)
     
subplot(3,2,4)
makeGammaFig(corrupt_pct,gamma_values) 

load('FA_output_Sim_1_k_1_comm_1_delta_5.mat','corrupt_pct','R_EM','R_REM','gamma_values'); 
 
subplot(3,2,5)
makeCongruenceFig(corrupt_pct,R_EM,R_REM)
     
subplot(3,2,6)
makeGammaFig(corrupt_pct,gamma_values) 

set(gcf,'Position',[50 50 600 800])

print('Fig3','-depsc','-r600','-painters')

function makeCongruenceFig(corrupt_pct,R_EM,R_REM)
       
    maj_pct = 1 - corrupt_pct;

    hold on
    errorbar(maj_pct, R_EM(:,1), R_EM(:,2), ...
        'LineWidth',1.5,'Color','b','LineStyle',':','Marker','d','MarkerSize',5)
    errorbar(maj_pct, R_REM(:,1), R_REM(:,2), ...
        'LineWidth',1.5,'Color','r','LineStyle','-','Marker','d','MarkerSize',5)
    set(gca,'YLim',[0.5 1])  
    ylabel({'RV';'Coefficient'})
    xlabel('Majority %')
    PrettyFig
    legend('EM','REM', 'Location','southeast','Orientation','horizontal')
    hold off
end

function makeGammaFig(corrupt_pct,gamma_values) 
    
    maj_pct = 1 - corrupt_pct;

    hold on
    errorbar(maj_pct, gamma_values(:,1), gamma_values(:,2),'LineWidth',1.5,'Color','m')
    plot(0:0.01:1, 0:0.01:1,'--','LineWidth',1,'Color','k')
    hold off
    ylabel({'Estimated';'Gamma'})
    xlabel('Majority %')
    set(gca,'YLim',[0.6 1])  
    PrettyFig
    legend off
end

end