clc
clear
rng(3)
betahat = randn(5,1);
k = randn(5,5);
sigma = k*k'+eye(5);
sigma = diag(sqrt(diag(sigma)))^-1*sigma*diag(sqrt(diag(sigma)))^-1;
sigma = sigma +diag(0:0.1:0.4);
Mbar = 1;
numPrePeriods = 4;
numPostPeriods = 1;
l_vec = 1;
alpha = 0.05;
hybrid_kappa = 0.005;
gridPoints = 1000;
grid_ub =  20;
grid_lb = -20;

bound = 'deviation from parallel trends';
method = 'C-LF';
[lb, ub] = createSensitivityResults_relativeMagnitudes(betahat, sigma, numPrePeriods, numPostPeriods, ...
    bound, method, Mbar, l_vec, alpha, gridPoints, nan, nan, hybrid_kappa);
[lb ub]
% %%
% CI_RM = computeConditionalCS_DeltaRM(betahat, sigma, numPrePeriods, numPostPeriods, l_vec, Mbar, alpha, hybrid_flag, hybrid_kappa, false, postPeriodMomentsOnly, gridPoints, grid_ub, grid_lb);
% %%
% thetaGrid = linspace(grid_lb, grid_ub, gridPoints);
% % resultsGrid = testOverThetaGrid(betahat, sigma,A_RM_s, d_RM, thetaGrid, numPrePeriods, alpha, @testInIdentifiedSet_LF_Hybrid, hybrid_list);
% % 
% % theta = -5.1;
% % a = testInIdentifiedSet_LF_Hybrid(betahat - basisVector(numPrePeriods + 1, length(betahat)) * theta, sigma, A_RM_s, d_RM, alpha, hybrid_list)
% 
% % theta = 3.8;
% % b = testInIdentifiedSet_LF_Hybrid(betahat - basisVector(numPrePeriods + 1, length(betahat)) * theta, sigma, A_RM_s, d_RM, alpha, hybrid_list)
% CI = testOverThetaGrid(betahat, sigma, A_RM_s, d_RM, thetaGrid, numPrePeriods, alpha, @testInIdentifiedSet_LF_Hybrid, hybrid_list);
% 
% %% why is it so slow
% 
% 
% 
% %%
% clear
% s = -2;
% max_positive = true;
% Mbar = 1;
% betahat = [-0.4277,-0.5794,0.9260, 0.0055, -0.6345]';
% sigma = eye(5);
% numPrePeriods = 4;
% numPostPeriods = 1;
% l_vec = 1;
% alpha = 0.05;
% hybrid_flag = 'LF';
% hybrid_kappa = 0.005;
% postPeriodMomentsOnly = true;
% gridPoints = 1000;
% grid_ub = 20;
% grid_lb = -20;
% 
% A_RM_s = create_A_RM(numPrePeriods, numPostPeriods, Mbar, s, max_positive);
% d_RM = create_d_RM(numPrePeriods, numPostPeriods);
% rowsForARP = 1:size(A_RM_s, 1);
% 
% % Create hybrid_list object
% hybrid_list = struct('hybrid_kappa', hybrid_kappa);
% if strcmp(hybrid_flag, 'LF')
%     % Compute LF CV and store it in hybrid_list
%     lf_cv = compute_least_favorable_cv([], A_RM_s * sigma * A_RM_s', hybrid_kappa);
%     hybrid_list.lf_cv = lf_cv;
% end
%     
% CI = ARP_computeCI(betahat, sigma, numPrePeriods, numPostPeriods, A_RM_s, ...
%     d_RM, l_vec, alpha, hybrid_flag, hybrid_list, false, grid_lb, grid_ub, gridPoints, rowsForARP);