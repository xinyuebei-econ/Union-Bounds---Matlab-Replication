function id_set = compute_IDset_DeltaRM(Mbar, trueBeta, l_vec, numPrePeriods, numPostPeriods)
    % This function computes the upper and lower bound of the identified set
    % given the event study coefficients, l_vec, and Mbar.
    %
    % To do so, we construct the identified set at each choice of s, +, -. We
    % then take the union of these intervals.
    %
    % Note: l_vec is assumed to be non-negative.
    %
    % Inputs:
    %   Mbar            = smoothness param of Delta^MB
    %   trueBeta        = vector of population event study coefficients
    %   l_vec           = vector l defining parameter of interest
    %   numPrePeriods   = number of pre-periods
    %   numPostPeriods  = number of post-periods
    %
    % Outputs:
    %   Struct with fields:
    %     id_ub = upper bound of ID set
    %     id_lb = lower bound of ID set

    % Construct identified sets for (+) at each value of s
    min_s = -(numPrePeriods - 1);
    id_bounds_plus = [];
    id_bounds_minus = [];
    for s = min_s:0
        id_bounds_plus_s = compute_IDset_DeltaRM_fixedS(s, Mbar, trueBeta, l_vec, numPrePeriods, numPostPeriods, true);
        id_bounds_minus_s = compute_IDset_DeltaRM_fixedS(s, Mbar, trueBeta, l_vec, numPrePeriods, numPostPeriods, false);
        id_bounds_plus = [id_bounds_plus; id_bounds_plus_s];
        id_bounds_minus = [id_bounds_minus; id_bounds_minus_s];
    end

    % Construct the identified set by taking the max of the upper bound and the min of the lower bound
    id_lb = min(min(id_bounds_plus.id_lb), min(id_bounds_minus.id_lb));
    id_ub = max(max(id_bounds_plus.id_ub), max(id_bounds_minus.id_ub));

    % Return identified set
    id_set = struct('id_lb', id_lb, 'id_ub', id_ub);
end


