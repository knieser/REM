function [R_EM, R_REM, gamma_values, eps_values, R_EM_minor, R_REM_minor] = FA_sim_main(sim,corrupt_pct,k_choice,delta)

% sim: (1) vary corruption (2) vary K
% can enter array into either corrupt_pct or k_choice

% e.g. FA_sim_main(1,0:0.1:0.4,1,0.05);


%%%%% Set Parameters %%%%%%

% # simulations to run to generate error bars;
n_sim = 25; 

% factor structure parameters;
sparse = 1; % if 1 will use Tucker-Koopman-Linn method to generate pop correlation matrices;
communality = 3; % 1 = low, 2 = wide, 3 = high, 4 = perfect;
m = 0; % controls similarity of lambdas, m=1 means the lambdas are the same;

% dimensions
p_0 = 10;
k_0 = 1;
n = 1000;

maj_pct = 1 - corrupt_pct;

%%%%%% Simulate population covariance/factor structures %%%%%%%

[sigma_01, sigma_02, lambda_01, lambda_02, ~] = FA_data_sim(p_0,k_0,sparse,communality,m);


%%%%%% Run simulation %%%%%%%

if sim == 1
    msg = 'Simulation 1: Varying Mix';
    disp(msg)
    
    % Initialize;
    R_EM = zeros(length(corrupt_pct), 2);
    R_REM = zeros(length(corrupt_pct), 2);
    
    R_EM_minor = zeros(length(corrupt_pct), 2);
    R_REM_minor = zeros(length(corrupt_pct), 2);
    
    gamma_values = zeros(length(corrupt_pct), 2);
    eps_values = zeros(length(corrupt_pct), 2);
    
    % Run through varying corrupt_pct
    for c = 1:length(corrupt_pct)
        msg = ['Working on corrupt pct = ', num2str(100*corrupt_pct(c)),'%'];
        disp(msg)
        
        R_classic = zeros(n_sim,1);
        R_robust = zeros(n_sim,1);
        R_classic_b = zeros(n_sim,1);
        R_robust_b = zeros(n_sim,1);
        gamma = zeros(n_sim,1);
        eps = zeros(n_sim,1);
        
        for l = 1:n_sim
            
            % Sample data from population;
            [X,lambda_target_1,lambda_target_2] = FA_sample_data(sigma_01,sigma_02,lambda_01,lambda_02,n,corrupt_pct(c));
            
            % Get estimates;
            [hlambda_1, ~, hlambda_2, ~, ~, gamma(l), eps(l)] = FA_estimates(X,k_choice,delta);
            
            % Compute and store congruence coefficient
            R_classic(l) = computeCongruence(hlambda_1, lambda_target_1); % first matrix gets procrustes rotation;
            R_robust(l) = computeCongruence(hlambda_2, lambda_target_1);
            
            R_classic_b(l) = computeCongruence(hlambda_1, lambda_target_2); % first matrix gets procrustes rotation;
            R_robust_b(l) = computeCongruence(hlambda_2, lambda_target_2);
          
        end
        
        % Aggregate data;
        R_EM(c,:) = aggregate(R_classic);
        R_REM(c,:) = aggregate(R_robust);
        
        R_EM_minor(c,:) = aggregate(R_classic_b);
        R_REM_minor(c,:) = aggregate(R_robust_b);
        
        gamma_values(c,:) = aggregate(gamma);
        eps_values(c,:) = aggregate(eps);
                     
    end 

    
    %%%%%%%%%%% Make Figures %%%%%%%%%%%%%%%%%%%

    % Plot congruence vs mix %;
    figure;
    hold on
    errorbar(maj_pct, R_EM(:,1), R_EM(:,2),'LineWidth',1.5,'Color','b')
    errorbar(maj_pct, R_REM(:,1), R_REM(:,2),'LineWidth',1.5,'Color','r')
    set(gca,'YLim',[0 1])  
    ylabel('Congruence Coeff')
    xlabel('Majority %')
    PrettyFig
    legend('Classic','Robust', 'Location','southeast','Orientation','vertical')
    hold off
    %print('Heterogeneity high communality m0','-dpng','-r300')
    
    
    % Plot congruence to minority vs mix %;
    figure;
    hold on
    errorbar(maj_pct, R_EM_minor(:,1), R_EM_minor(:,2),'LineWidth',1.5,'Color','b')
    errorbar(maj_pct, R_REM_minor(:,1), R_REM_minor(:,2),'LineWidth',1.5,'Color','r')
    set(gca,'YLim',[0 1])  
    ylabel('Congruence Coeff')
    xlabel('Majority %')
    PrettyFig
    legend('Classic','Robust', 'Location','southeast','Orientation','vertical')
    hold off
    %print('Heterogeneity high communality m0','-dpng','-r300')
    
    
    % Plot estimated gamma vs mix %;
    figure
    hold on
    errorbar(maj_pct, gamma_values(:,1), gamma_values(:,2),'LineWidth',1.5,'Color','m')
    plot(0:0.01:1, 0:0.01:1,'--','LineWidth',1,'Color','k')
    hold off
    ylabel('Estimated Gamma')
    xlabel('Majority %')
    set(gca,'YLim',[0.6 1])  
    PrettyFig
    legend('off')
    %print('Gamma ep16 n2000','-dpng','-r300')
    
elseif sim == 2
    msg = 'Simulation 2: Varying K';
    disp(msg)
    
    % Initialize;
    R_EM = zeros(length(k_choice), 2);
    R_REM = zeros(length(k_choice), 2);
    
    R_EM_minor = zeros(length(k_choice), 2);
    R_REM_minor = zeros(length(k_choice), 2);
    
    gamma_values = zeros(length(k_choice), 2);
    eps_values = zeros(length(k_choice), 2);
         
    % Run through varying k_choice
    for k = 1:length(k_choice)
        msg = ['Working on k = ', num2str(k_choice(k))];
        disp(msg)
        
        R_classic = zeros(n_sim,1);
        R_robust = zeros(n_sim,1);
        R_classic_b = zeros(n_sim,1);
        R_robust_b = zeros(n_sim,1);
        gamma = zeros(n_sim,1);
        eps = zeros(n_sim,1);

        for l = 1:n_sim
        
            % Sample data from population;
            [X,lambda_target_1,lambda_target_2] = FA_sample_data(sigma_01, sigma_02, lambda_01, lambda_02,n,corrupt_pct);
            
            % Get estimates;
            [hlambda_1, ~, hlambda_2, ~,  ~, gamma(l), eps(l)] = FA_estimates(X,k_choice(k),delta);
            
            % Compute and store congruence coefficient
            R_classic(l) = computeCongruence(hlambda_1, lambda_target_1); % first matrix gets procrustes rotation;
            R_robust(l) = computeCongruence(hlambda_2, lambda_target_1);
            
            R_classic_b(l) = computeCongruence(hlambda_1, lambda_target_2); % first matrix gets procrustes rotation;
            R_robust_b(l) = computeCongruence(hlambda_2, lambda_target_2);
        
        end

        % Aggregate data;
        R_EM(k,:) = aggregate(R_classic);
        R_REM(k,:) = aggregate(R_robust);
        
        R_EM_minor(k,:) = aggregate(R_classic_b);
        R_REM_minor(k,:) = aggregate(R_robust_b);
        
        gamma_values(k,:) = aggregate(gamma);
        eps_values(k,:) = aggregate(eps);
        

    end 
    
    %%%%%%%%%%% Make Figures %%%%%%%%%%%%%%%%%%%
    
    % Plot congruence vs mix %;
    figure;
    hold on
    errorbar(k_choice, R_EM(:,1), R_EM(:,2),'LineWidth',1.5,'Color','b')
    errorbar(k_choice, R_REM(:,1), R_REM(:,2),'LineWidth',1.5,'Color','r')
    set(gca,'YLim',[0 1])  
    ylabel('Congruence Coeff')
    xlabel('K Choice')
    PrettyFig
    legend('Classic','Robust', 'Location','southeast','Orientation','vertical')
    hold off
    %print('Heterogeneity high communality m0','-dpng','-r300')
    
    
    % Plot estimated gamma vs mix %;
    figure
    hold on
    errorbar(k_choice, gamma_values(:,1), gamma_values(:,2),'LineWidth',1.5,'Color','m')
    hold off
    ylabel('Estimated Gamma')
    xlabel('K Choice')
    set(gca,'YLim',[0 1])  
    PrettyFig
    legend('off')
    %print('Gamma ep16 n2000','-dpng','-r300')
end
 
save FA_sim_5_14_20;

end



