% simulation
clc
clear

addpath('./data')
addpath('./aux_RR23')
addpath('./aux_YKHS')
addpath('./aux_ModiCon')

% see calibrateDGP.m for details
Type = 1;
DGP = 6;

M = 1;
alpha  = 0.1;
alphac = alpha*4/5;
[delta0, deltaSigma, T_pre, Al, Au, index] = calibrateDGP(DGP, Type, M);
filename = strcat('./results/simulation_DGP',num2str(DGP),'_Type',num2str(Type));

n = 5000;
B = 2500;
S = 1000;
Blarge = B*10;
eta    = 0.001;

lambdaSigma = [Al; Au]*deltaSigma*[Al; Au]';
lambdasigma = sqrt(diag(lambdaSigma));
kk = size(Al,1);
sigma_l  = lambdasigma(1:kk);
sigma_u  = lambdasigma(kk+1:end);
tol    = 10^-3;
tol_r  = 10^-3;

l_vec = 1;
hybrid_kappa = alpha/10;
gridPoints = 2000;
bound = 'deviation from parallel trends';
method = 'C-LF';

num = feature('numCores');
parpool(num)
parfor s = 1 : S
    s
    tic
    rng(s,'twister')
    data = mvnrnd(delta0, deltaSigma*n, n);    
    deltahat  = mean(data)';
 
    lambdahat_l = Al*deltahat;
    lambdahat_u = Au*deltahat; 
    
    [CI_h(s,:), CI_c(s,:), CI_p(s,:), delta1(:,s), c(s)] = CIhybrid(deltahat, deltaSigma, Al, Au, alpha, alphac, eta, B, Blarge, tol, tol_r, index);
    CI_sim(s,:)  = [min(lambdahat_l - norminv(1-alpha/2)*sigma_l) max(lambdahat_u + norminv(1-alpha/2)*sigma_u)];
    CI_OS(s,:)   = [min(lambdahat_l - norminv(1-alpha)*sigma_l) max(lambdahat_u + norminv(1-alpha)*sigma_u)];
    CI_pt(s,:)   = [min(lambdahat_l) max(lambdahat_u)];
    grid_ub   = min(lambdahat_l - 5*sigma_l);
    grid_lb   = max(lambdahat_u + 5*sigma_u);    
    [lb, ub]  = createSensitivityResults_relativeMagnitudes(deltahat, deltaSigma, T_pre, 1, ...
                                    bound, method, M, l_vec, alpha, gridPoints, grid_ub, grid_lb, hybrid_kappa);
    CI_RR(s,:) = [lb, ub]; 
    
    % Bootstrap & adjusted bootstrap
    m = ceil(n/log(log(n)));
    [CI_YKHS_set(s,:), CI_YKHS_pt(s,:)] = CI_boot_pt(deltahat, data, deltaSigma, Al, Au, alpha, n, m, B);
    time(s) = toc;
end
save(filename)
