function A = create_A_RM(numPrePeriods, numPostPeriods, Mbar, s, max_positive, dropZero)
    % This function creates a matrix for the linear constraints that
    % delta is in Delta^RM_{s,(.)}(Mbar), where (.) is + if max_positive = T and (-) if max_positive = F.
    %
    % Inputs:
    %   numPrePeriods  = number of pre-periods.
    %   numPostPeriods = number of post-periods.
    %   Mbar           = smoothness parameter of Delta^MB. Default is 1.
    %   s              = parameter determining the position of the first difference. Default is 0.
    %   max_positive   = flag indicating whether to use positive moments or negative moments. Default is true.
    %   dropZero       = flag to drop the constraint corresponding to t = 0. Default is true.

    if nargin < 3
        Mbar = 1;
    end
    if nargin < 4
        s = 0;
    end
    if nargin < 5
        max_positive = true;
    end
    if nargin < 6
        dropZero = true;
    end

    % First construct matrix Atilde that takes first differences
    Atilde = zeros(numPrePeriods + numPostPeriods, numPrePeriods + numPostPeriods + 1);
    for r = 1 : (numPrePeriods + numPostPeriods)
        Atilde(r, r : (r + 1)) = [-1, 1];
    end

    % Create a vector to extract the max first difference, which corresponds with the first difference for period s,
    % or minus this if max_positive == false
    v_max_dif = zeros(1, numPrePeriods + numPostPeriods + 1);
    v_max_dif((numPrePeriods + s ) : (numPrePeriods + s + 1)) = [-1, 1];

    if ~max_positive
        v_max_dif = -v_max_dif;
    end

    % The bounds for the first difference starting with period t are 1 * v_max_dif if t <= 0
    % and M * v_max_dif if t > 0
    A_UB = [repmat(v_max_dif, numPrePeriods, 1); repmat(Mbar * v_max_dif, numPostPeriods, 1)];

    % Construct A that imposes |Atilde * delta | <= A_UB * delta
    A = [Atilde - A_UB; -Atilde - A_UB];

    % Remove all-zero rows of the matrix Atilde
    zerorows = sum(A .^ 2, 2) <= 1e-10;
    A = A(~zerorows, :);

    % Remove the period corresponding to t = 0
    if dropZero        
        A(:, numPrePeriods + 1) = [];
    end
end


