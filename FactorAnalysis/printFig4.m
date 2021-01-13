function printFig4(delta)

close all;

set(0,'defaulttextinterpreter','latex')

load(['FA_output_Sim_1_comm_3_delta_', num2str(100*delta),'.mat'],'corrupt_pct','R_EM','R_REM','gamma_values'); 
 
subplot(3,2,1)
makeCongruenceFig(corrupt_pct,R_EM,R_REM)
ylabel({'$$\textbf{High}$$';'';'';'Congruence'});
     
subplot(3,2,2)
makeGammaFig(corrupt_pct,gamma_values) 
ylabel('Est. Gamma');

load(['FA_output_Sim_1_comm_2_delta_', num2str(100*delta),'.mat'],'corrupt_pct','R_EM','R_REM','gamma_values'); 
 
subplot(3,2,3)
makeCongruenceFig(corrupt_pct,R_EM,R_REM)
ylabel({'$$\textbf{Wide}$$';'';'';'Congruence'});
     
subplot(3,2,4)
makeGammaFig(corrupt_pct,gamma_values) 
ylabel('Est. Gamma');

load(['FA_output_Sim_1_comm_1_delta_', num2str(100*delta),'.mat'],'corrupt_pct','R_EM','R_REM','gamma_values'); 
 
subplot(3,2,5)
makeCongruenceFig(corrupt_pct,R_EM,R_REM)
ylabel({'$$\textbf{Low}$$';'';'';'Congruence'}); 
xlabel('Majority, proportion');

     
subplot(3,2,6)
makeGammaFig(corrupt_pct,gamma_values) 
ylabel('Est. Gamma'); 
xlabel('Majority, proportion');


set(gcf,'Position',[50 50 800 1000])


print('Fig4','-depsc','-r600','-painters')

function makeCongruenceFig(corrupt_pct,R_EM,R_REM)
    set(0,'defaulttextinterpreter','latex')
   
    maj_pct = 1 - corrupt_pct;

    hold on
    errorbar(maj_pct, R_EM(:,1), R_EM(:,2), ...
        'LineWidth',1.5,'Color','b','LineStyle','--','Marker','d','MarkerSize',5)
    errorbar(maj_pct, R_REM(:,1), R_REM(:,2), ...
        'LineWidth',1.5,'Color','r','LineStyle','-','Marker','d','MarkerSize',5)
    set(gca,'YLim',[0.8 1])  
    set(gca,'XLim',[0.6 1]) 
    xticks(0.6:0.1:1)
    yticks(0.8:0.05:1)
    PrettyFig
    legend('EM','REM', 'Location','southeast','Orientation','vertical','FontSize',9)
    hold off
end

function makeGammaFig(corrupt_pct,gamma_values)  
    set(0,'defaulttextinterpreter','latex')

    maj_pct = 1 - corrupt_pct;

    hold on
    errorbar(maj_pct, gamma_values(:,1), gamma_values(:,2),'LineWidth',1.5,'Color','m')
    plot(0:0.01:1, 0:0.01:1,'--','LineWidth',1,'Color','k')
    hold off
    set(gca,'YLim',[0.6 1])  
    set(gca,'XLim',[0.6 1])
    xticks(0.6:0.1:1)
    PrettyFig
    legend off  
end

end