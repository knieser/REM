%%%%%%%%%% Main Script to Run All Simulations %%%%%%%%%%%%

% Select delta
delta = 0.05;

% Run simulation for Example 1
GMM_sim_main(2,2,delta);

% Check AIC/BIC for Example 1
[AIC_sim_2_EM, BIC_sim_2_EM, AIC_sim_2_REM, BIC_sim_2_REM] = selectModel(2,delta);

% Run simulation for Example 2
GMM_sim_main(3,2,delta);

% Check AIC/BIC for Example 2
[AIC_sim_3_EM, BIC_sim_3_EM, AIC_sim_3_REM, BIC_sim_3_REM]= selectModel(3,delta);

save('model_selection_output.mat');

% Run simulations with different K (Fig 3)
GMM_sim_main(1,1,delta);
GMM_sim_main(1,2,delta);
GMM_sim_main(1,3,delta);

% Extra simulations
GMM_sim_main(4,2,delta);
GMM_sim_main(5,2,delta);
GMM_sim_main(6,2,delta);

% Print figures;
printFig1(delta); close all;
printFig2(delta); close all;
printFig3(delta); close all;
printFig6(delta); close all;
printFig7(delta); close all;
