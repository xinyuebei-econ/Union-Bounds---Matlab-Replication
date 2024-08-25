function CI = computeConditionalCS_DeltaRM_fixedS(s, max_positive, Mbar, betahat, sigma, numPrePeriods,numPostPeriods, l_vec, alpha, hybrid_flag, hybrid_kappa, postPeriodMomentsOnly, gridPoints, grid_ub, grid_lb)
% This function computes the ARP CI that includes nuisance parameters
% for Delta^{RM}(Mbar) for a fixed s and (+),(-). This function uses ARP_computeCI for all
% of its computations. It is used as a helper function in computeConditionalCS_DeltaRM below.
% Check that hybrid_flag equals LF or ARP
if ~strcmp(hybrid_flag, 'LF') && ~strcmp(hybrid_flag, 'ARP')
    error('hybrid_flag must equal ''ARP'' or ''LF''');
end

% Create hybrid_list object
hybrid_list = struct('hybrid_kappa', hybrid_kappa);

% Create matrix A_RM_s, and vector d_RM
A_RM_s = create_A_RM(numPrePeriods, numPostPeriods, Mbar, s, max_positive);
d_RM = create_d_RM(numPrePeriods, numPostPeriods);

% If only use post period moments, construct indices for the post period moments only.
if postPeriodMomentsOnly && (numPostPeriods > 1)
    postPeriodIndices = (numPrePeriods + 1):size(A_RM_s, 2);
    postPeriodRows = find(sum(A_RM_s(:, postPeriodIndices) ~= 0, 2) > 0);
    rowsForARP = postPeriodRows;
else
    rowsForARP = 1:size(A_RM_s, 1);
end

% if there is only one post-period, we use the no-nuisance parameter functions
if numPostPeriods == 1
    if strcmp(hybrid_flag, 'LF')
        % Compute LF CV and store it in hybrid_list
        lf_cv = compute_least_favorable_cv([], A_RM_s * sigma * A_RM_s', hybrid_kappa);
        hybrid_list.lf_cv = lf_cv;
    end
    
    % Compute confidence set
    CI = APR_computeCI_NoNuis(betahat, sigma, A_RM_s, d_RM, numPrePeriods, numPostPeriods, l_vec, alpha, false, hybrid_flag, hybrid_list, grid_ub, grid_lb, gridPoints);
else
    % Compute ARP CI for l'beta using Delta^RM
    warning('only for numPost = 1')
    CI = ARP_computeCI(betahat, sigma, numPrePeriods, numPostPeriods, A_RM_s, d_RM, l_vec, alpha, hybrid_flag, hybrid_list, false, grid_lb, grid_ub, gridPoints, rowsForARP);
        
end
end

