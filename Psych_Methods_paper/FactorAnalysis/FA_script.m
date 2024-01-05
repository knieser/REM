
% Select delta
delta = 0.05;


% Run simulations changing corruption at each communality level
tic
FA_sim_main(1,3,delta);
toc
tic
FA_sim_main(1,2,delta);
toc
tic
FA_sim_main(1,1,delta);
toc
% Run simulations for changing K (with high communality)
tic
FA_sim_main(2,3,delta);
toc
% Print figures;
printFig4(delta);
printFig5(delta);
