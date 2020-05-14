function [est_gmm,est_gmm_2,est_weights,est_gamma,opt_eps] = LPA_sim_main(X,k_choice,global_iter,delta,step,seed)

rng(seed)
    
% Set # iterations for epsilon search;
maxiter = step+1;
    

%%%%%%% Initialize %%%%%%%%%%%%

est_mu = cell(length(k_choice),1);
est_sigma = cell(length(k_choice),1);
est_mix = cell(length(k_choice),1);


%%%%%%%% Estimation %%%%%%%%%%%%%

for k = 1:length(k_choice)

     % Get estimates from EM; 
     
    disp('Starting global opt for EM...')
       
    [est_mu{k}, est_sigma{k}, est_mix{k}, ~, ~, ~, ~, ind_llh] ...
        = globalRobustEMAlgMixture(X,k_choice(k),global_iter,0,0);  
    
    disp('...done')
    
    
    % Get estimates from REM;
    chk = zeros(1,maxiter);
    est_gamma = zeros(1,maxiter);
    
    eps_max = quantile(exp(ind_llh), 0.50); % Median likelihood value;
    eps_range = 0:eps_max/step:eps_max;
    
    
    % Get initial seeds;
    GMModel = fitgmdist(X',k_choice(k),'Options',statset('MaxIter',500));
    intl_sigma = GMModel.Sigma;
    intl_mix = GMModel.ComponentProportion;
       
    % Global optimization to get starting mu;
    eps_start = quantile(exp(ind_llh),0.05);
    disp('Starting global opt for REM...')
    [intl_mu] = globalRobustEMAlgMixture(X, k_choice(k), global_iter, eps_start,1);
    disp('...done')
 
    for iter = 1:maxiter
    
        % Get REM estimates;
        [mu, sigma, mix, ~, est_gamma(iter), est_weights] ...
            = RobustEMAlgMixture(X,k_choice(k),eps_range(iter), intl_mu, intl_sigma, intl_mix);
        
        % Calculate prob. that weights <= 1/2;
        chk(iter) = checkEps(mu,sigma,mix,est_gamma(iter),eps_range(iter));
                
        if chk(iter) > delta 
            break;
        end
    end
    
    % Truncate gamma and chk;
    chk = chk(1:iter);
    est_gamma = est_gamma(1:iter);
               
    opt_eps = eps_range(iter);
    
end 

est_gmm = gmdistribution(est_mu{k}', est_sigma{k}, est_mix{k});
est_gmm_2 = gmdistribution(mu', sigma, mix);


    
end