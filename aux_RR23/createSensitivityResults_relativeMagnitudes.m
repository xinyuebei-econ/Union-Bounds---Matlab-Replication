function [lb, ub] = createSensitivityResults_relativeMagnitudes(betahat, sigma, numPrePeriods, numPostPeriods, bound, method, Mbarvec, l_vec, alpha, gridPoints, grid_ub, grid_lb, hybrid_kappa)
% rewrote by ChatGPT 3.5
% If Mbarvec is null, construct default Mbarvec to be 10 values on [0,2].
    if isempty(Mbarvec)
        Mbarvec = linspace(0, 2, 10);
    end

    % Use method to specify the choice of confidence set
    if strcmp(method, 'C-LF')
        hybrid_flag = 'LF';
        method_named = 'C-LF';
    elseif strcmp(method, 'Conditional')
        hybrid_flag = 'ARP';
        method_named = 'Conditional';
    else
        error('method must be either NULL, Conditional or C-LF.');
    end

    % If bound = "parallel trends violation", we select Delta^{RM} and its variants.
    if strcmp(bound, 'deviation from parallel trends')
        % use monotonicity direction and biasDirection to select the choice of Delta
        % Mbarvec is a scalar here
        temp = computeConditionalCS_DeltaRM(betahat, sigma, numPrePeriods, numPostPeriods, l_vec, Mbarvec, ...
            alpha, hybrid_flag, hybrid_kappa, false, true, gridPoints, grid_ub, grid_lb);
        lb = min(temp(temp(:,2) == 1,1));
        ub = max(temp(temp(:,2) == 1,1));
        
  
    else % if bound = "deviation from linear trend", we select Delta^{SDRM} and its variants.
        error('The code only includes RM part')      
    end
end


