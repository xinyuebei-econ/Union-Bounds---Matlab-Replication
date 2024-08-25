function CI = APR_computeCI_NoNuis(betahat, sigma, A, d, numPrePeriods, numPostPeriods, l_vec, alpha, returnLength, hybrid_flag, hybrid_list, grid_ub, grid_lb, gridPoints)
    % This function computes the APR confidence interval for Delta^SD(M) for a given event study.
    % It takes the following inputs:
    %   betahat = vector of estimated event study coefficients
    %   sigma   = covariance matrix of estimated event study coefficients
    %   A       = matrix defining the set Delta
    %   d       = vector defining the set Delta
    %   numPrePeriods = number of pre-periods
    %   numPostPeriods = number of post-periods
    %   M = tuning parameter of Delta^SD(M), default M = 0
    %   postPeriod = post-period of interest
    %   alpha = size of CI, default 0.05.
    
    % Construct grid of tau values to test and test which values of tau lie in ID set.
    thetaGrid = linspace(grid_lb, grid_ub, gridPoints);
    
    if strcmp(hybrid_flag, 'ARP')
        resultsGrid = testOverThetaGrid(betahat, sigma, A, d, thetaGrid, numPrePeriods, alpha, @testInIdentifiedSet);
    elseif strcmp(hybrid_flag, 'FLCI')
        warning('not FLCI')
        resultsGrid = testOverThetaGrid(betahat, sigma, A, d, thetaGrid, numPrePeriods, alpha, @testInIdentifiedSet_FLCI_Hybrid, hybrid_list);
                      
    elseif strcmp(hybrid_flag, 'LF')
        resultsGrid = testOverThetaGrid(betahat, sigma, A, d, thetaGrid, numPrePeriods, alpha, @testInIdentifiedSet_LF_Hybrid, hybrid_list);
    else
        error('hybrid_flag must equal ''APR'' or ''FLCI'' or ''LF''');
    end
    resultsGrid(:,2) = 1-resultsGrid(:,2);
    if resultsGrid(1,2) == 1 || resultsGrid(end,2) == 1
        warning('CI is open at one of the endpoints; CI length may not be accurate');
    end
    
    % Compute length, else return grid
    if returnLength
        gridLength = 0.5 * ( [0, diff(thetaGrid)] + [diff(thetaGrid), 0] );
        CI = sum(resultsGrid(:, 2) .* gridLength);
    else
        CI = table(resultsGrid(:, 1), resultsGrid(:, 2), 'VariableNames', {'grid', 'accept'});
    end
end


