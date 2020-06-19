function printFig5()

clear; close all;

set(0,'defaulttextinterpreter','latex')

load('FA_output_Sim_2_corrupt_30pct_delta_5.mat','k_0','k_choice','R_EM','R_REM'); 

hold on
errorbar(k_choice, R_EM(:,1), R_EM(:,2),'LineWidth',1.5,'Color','b','LineStyle','--','Marker','d','MarkerSize',5)
errorbar(k_choice, R_REM(:,1), R_REM(:,2),'LineWidth',1.5,'Color','r','LineStyle','-','Marker','d','MarkerSize',5)
xline(k_0,':','label',{'True Number';'of Factors'}, ...
    'FontSize', 14, 'LabelVerticalAlignment','middle','LabelOrientation','horizontal');
xticks(1:1:6)
yticks(0.8:0.05:1)
set(gca,'YLim',[0.8 1])  
set(gca,'XLim',[1,6])
ylabel('Congruence')
xlabel('Selected Number of Factors')
PrettyFig
legend('EM','REM', 'Location','southeast','Orientation','vertical')
title('$$\textbf{Varying Number of Factors}$$')
hold off

set(gcf,'Position',[50 50 800 600])

print('Fig5','-depsc','-r600','-painters')

%{
y0 = cumsum(mean(lambda_01.^2))';
y1 = mean(avg_EM_h2)';
y2 = mean(avg_REM_h2)';

bar(k_choice',[y0,y1,y2]);
PrettyFig
legend('Truth','EM','REM','location','southoutside','orientation','horizontal','FontSize',10)
xlabel('Factors');
ylabel({'Cumulative';'% Variance';'Explained'}); ylim([0,1]); yticks(0:0.2:1);
%}
end