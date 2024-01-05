function [X, true_gmm, grp_flag] = GMM_data_sim(sim,p,k,n,mix)
%{
This is the function generates the data.
    
INPUT:
    sim_num: 
        1 = No corruption,
        2 = Scattered minority,
        3 = Scattered minority w cross
    p: number of observed variables
    k: number of latent groups
    n: number of observations
    mix: (1 x k) vector of mixture probabilities 
    
OUTPUT:
    X: (p x n) array of simulated data
    true_gmm: ground truth distribution
    grp_flag: (n x 1) vector of group indicators
%}

% Set parameter values based on simulation number
if sim == 1
    msg = 'Simulation 1: No Corruption';
    disp(msg)
    
    mu = [2,7,6; 7,7,3];
    sigma = reshape(1/8*[2, -1, 2, 1, 4, 0; -1, 2, 1, 2, 0, 4],p,p,[]);
    skew = 0.5; 
    
elseif sim == 2
    msg = 'Simulation 2: Scattered Minority';
    disp(msg)

    mu = [2.5,7,7; 7,2.5,7]; 
    sigma = reshape(1/2*[2, -1, 2, 1, 20, 0; -1, 2, 1, 2, 0, 20],p,p,[]);
    skew = 0.5;

elseif sim == 3
    msg = 'Simulation 3: Scattered Minority with Cross';
    disp(msg)
    
    mu =  [5,5,5; 5,5,6]; 
    sigma = reshape(1/2*[2, -1.6, 2, 1.6, 10, 0; -1.6, 2, 1.6, 2, 0, 10],p,p,[]);
    skew = 0.5; 
    
elseif sim == 4
    msg = 'Simulation 4: Scattered Minority, Large Within-Group Variance';
    disp(msg)
    
    mu = [3,6,7; 7,3,7];
    sigma = reshape(1/2*[2, 1, 6, 1, 20, 0; 1, 2, 1, 10, 0, 20],p,p,[]);
    skew = 0.5;

elseif sim == 5
    msg = 'Simulation 5: Scattered Minority, Large Within-Group Variance';
    disp(msg)
    
    mu = [3,6,7; 7,3,7];
    sigma = reshape(1/2*[2, 1, 6, 1, 20, 0; 1, 2, 1, 6, 0, 20],p,p,[]);
    skew = 0.5;
    
elseif sim == 6
    msg = 'Simulation 6: Noise within a Group';
    disp(msg)
    
    mu = [3,7,5; 7,3,2];
    sigma = reshape(1/2*[4, 1, 4, 1, 0.2, 0; 1, 4, 1, 4, 0, 0.2],p,p,[]);
    skew = 0.5;
end

% Make ground truth distribution
true_gmm = gmdistribution(mu',sigma,mix);

%%%%%%%% Simulate Data %%%%%%%%%%

% Pre-allocate;
R = cell(k,1);
X = zeros(p,n);

% Calculate Cholesky for each group's covariance matrix
for j = 1:k
    R{j} = chol(sigma(:,:,j));
end

% Choose random number for each data point;
U = rand(n,1);

% Flag to keep track of group;
grp_flag = zeros(n,1);

% Use random number to determine which distribution they came from;
for i = 1:n
    if U(i) < mix(1)
        X(:,i) = mu(:,1) + R{1}'*pearsrnd(0,1,skew,3,p,1);
        grp_flag(i) = 1;
    elseif U(i) < mix(1) + mix(2)
        X(:,i) = mu(:,2) + R{2}'*pearsrnd(0,1,skew,3,p,1);
        grp_flag(i) = 2;
    else
        if sim == 1
            X(:,i) = mu(:,3) + R{3}'*pearsrnd(0,1,skew,3,p,1);
        elseif sim == 2
            X(:,i) = 10*betarnd(1/2,1/3,p,1);
        elseif sim == 3
            X(:,i) = 10*betarnd(1/2,1/3,p,1);
        elseif sim == 4
            X(:,i) = 10*betarnd(1/2,1/3,p,1);
        elseif sim == 5
            X(:,i) = 10*betarnd(1/2,1/3,p,1);
        elseif sim == 6
            X(:,i) = mu(:,3) + R{3}'*pearsrnd(0,1,skew,3,p,1);
        end        
        grp_flag(i) = 3;
    end
end

end