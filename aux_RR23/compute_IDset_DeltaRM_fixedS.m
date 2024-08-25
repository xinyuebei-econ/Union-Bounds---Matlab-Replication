function id_bounds = compute_IDset_DeltaRM_fixedS(s, Mbar, max_positive, trueBeta, l_vec, numPrePeriods, numPostPeriods)
    % This helper function computes the upper and lower bound of the identified set
    % given the event study coefficients, l_vec, and Mbar. It computes the identified
    % set for a user-specified choice of s and (+), (-). This is used by
    % the function compute_IDset_DeltaRM below.

    % Create objective function: Wish to min/max l'delta_post
    fDelta = [zeros(1, numPrePeriods), l_vec];

    % Create A_RM, d_RM for this choice of s, max_positive
    A_RM_s = create_A_RM(numPrePeriods, numPostPeriods, Mbar, s, max_positive);
    d_RM = create_d_RM(numPrePeriods, numPostPeriods);

    % Create vector for direction of inequalities associated with RM
    dir_RM = repmat('<=', 1, length(d_RM));

    % Add equality constraint for pre-period coefficients
    prePeriodEqualityMat = [diag(ones(1, numPrePeriods)), zeros(numPrePeriods, numPostPeriods)];
    A_RM_s = [A_RM_s; prePeriodEqualityMat];
    d_RM = [d_RM; trueBeta(1:numPrePeriods)];
    dir_RM = [dir_RM, repmat('==', 1, size(prePeriodEqualityMat, 1))];

    % Specify variables between (-inf, inf)
    bounds.lower.ind = 1:(numPrePeriods + numPostPeriods);
    bounds.lower.val = -Inf(1, numPrePeriods + numPostPeriods);
    bounds.upper.ind = 1:(numPrePeriods + numPostPeriods);
    bounds.upper.val = Inf(1, numPrePeriods + numPostPeriods);

    % Create and solve for max
    save('lp','A_RM_s','d_RM')
    stop
    results_max = linprog(fDelta, A_RM_s, d_RM, [], [], bounds.lower.val, bounds.upper.val);

    % Create and solve for min
    results_min = linprog(-fDelta, A_RM_s, d_RM, [], [], bounds.lower.val, bounds.upper.val);

    if results_max.status ~= 1 && results_min.status ~= 1
        % If the solver does not return a solution, we just return the l_vec'trueBeta.
        id_ub = l_vec * trueBeta((numPrePeriods+1):(numPrePeriods+numPostPeriods))';
        id_lb = l_vec * trueBeta((numPrePeriods+1):(numPrePeriods+numPostPeriods))';
    else
        % Construct upper/lower bound of identified set
        id_ub = l_vec * trueBeta((numPrePeriods+1):(numPrePeriods+numPostPeriods))' - results_min.fval;
        id_lb = l_vec * trueBeta((numPrePeriods+1):(numPrePeriods+numPostPeriods))' - results_max.fval;
    end

    id_bounds = struct('id_lb', id_lb, 'id_ub', id_ub);
end


