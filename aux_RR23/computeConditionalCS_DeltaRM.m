function CI_RM = computeConditionalCS_DeltaRM(betahat, sigma, numPrePeriods, numPostPeriods, l_vec, Mbar, alpha, hybrid_flag, hybrid_kappa, returnLength, postPeriodMomentsOnly, gridPoints, grid_ub, grid_lb)
    % This function computes the ARP CI that includes nuisance parameters
    % for Delta^{RM}(Mbar). This function uses ARP_computeCI for all
    % of its computations.
    %
    % Inputs:
    %   betahat                 = vector of estimated event study coefficients
    %   sigma                   = covariance matrix of estimated event study coefficients
    %   numPrePeriods           = number of pre-periods
    %   numPostPeriods          = number of post-periods
    %   l_vec                   = vector that defines parameter of interest
    %   Mbar                    = tuning parameter for Delta^RM(Mbar), default Mbar = 0.
    %   alpha                   = desired size of CI, default alpha = 0.05
    %   hybrid_flag             = flag for hybrid, default = "LF". Must be either "LF" or "ARP"
    %   hybrid_kappa            = desired size of first-stage hybrid test, default = NULL
    %   returnLength            = returns length of CI only. Otherwise, returns matrix with grid in col 1 and test result in col 2.
    %   postPeriodMomentsOnly   = exclude moments for delta^MB that only include pre-period coefs
    %   gridPoints              = number of gridpoints to test over, default = 1000
    %   grid_ub                 = upper bound of the grid, default = NA
    %   grid_lb                 = lower bound of the grid, default = NA
    %
    % Outputs:
    %   data_frame containing upper and lower bounds of CI.

    % Create minimal s index for looping.
    min_s = -(numPrePeriods - 1);
    s_indices = min_s:0;

    % If grid_ub, grid_lb is not specified, we set these bounds to be equal to the id set under parallel trends
    % {0} +- 20*sdTheta (i.e. [-20*sdTheta, 20*sdTheta].
    sdTheta = sqrt(l_vec' * sigma((numPrePeriods+1):(numPrePeriods+numPostPeriods), (numPrePeriods+1):(numPrePeriods+numPostPeriods)) * l_vec);
    if isnan(grid_ub)
        grid_ub = 10 * sdTheta;
    end
    if isnan(grid_lb)
        grid_lb = -10 * sdTheta;
    end

    % Loop over s values for (+), (-), left join the resulting CIs based on the grid
    CIs_RM_plus_allS = zeros(gridPoints, length(s_indices));
    CIs_RM_minus_allS = zeros(gridPoints, length(s_indices));
    for s_i = 1:length(s_indices)
        % Compute CI for s, (+) and bind it to all CI's for (+)
        CI_s_plus = computeConditionalCS_DeltaRM_fixedS(s_indices(s_i), true, Mbar, betahat, sigma, numPrePeriods, numPostPeriods, l_vec, alpha, hybrid_flag, hybrid_kappa, postPeriodMomentsOnly, gridPoints, grid_ub, grid_lb);
                    
        CIs_RM_plus_allS(:, s_i) = CI_s_plus.accept;

        % Compute CI for s, (-) and bind it to all CI's for (-)
        CI_s_minus = computeConditionalCS_DeltaRM_fixedS(s_indices(s_i), false, Mbar, betahat, sigma, numPrePeriods, numPostPeriods, l_vec, alpha, hybrid_flag, hybrid_kappa, postPeriodMomentsOnly, gridPoints, grid_ub, grid_lb);
        CIs_RM_minus_allS(:, s_i) = CI_s_minus.accept;
    end

    CIs_RM_plus_maxS = max(CIs_RM_plus_allS, [], 2);
    CIs_RM_minus_maxS = max(CIs_RM_minus_allS, [], 2);

    % Take the max between (+), (-) and Construct grid containing theta points and whether any CI accepted
    grid = linspace(grid_lb, grid_ub, gridPoints)';
    accept = max(CIs_RM_plus_maxS, CIs_RM_minus_maxS);

    % Compute length if returnLength == true, else return grid
    if returnLength
        gridLength = 0.5 * ( [0; diff(grid)] + [diff(grid); 0] );
        CI_RM = sum(accept .* gridLength);
    else
        CI_RM = [grid, accept];
    end
end


