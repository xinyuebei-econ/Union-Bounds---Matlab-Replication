function d = create_d_RM(numPrePeriods, numPostPeriods, dropZero)
    % This function creates a vector for the linear constraints that
    % delta is in Delta^RM_{s,(.)}(Mbar), where (.) is + if max_positive = T and - if max_positive = F.
    % It implements this using the general characterization of d, NOT the sharp
    % characterization of the identified set.
    %
    % Inputs:
    %   numPrePeriods  = number of pre-periods.
    %   numPostPeriods = number of post-periods.
    %   dropZero       = flag to drop zero elements from d. Default is true.

    if nargin < 3
        dropZero = true;
    end

    A_RM = create_A_RM(numPrePeriods, numPostPeriods, 0, 0, dropZero); % d doesn't depend on Mbar or s; we just use this to get the dims right
    d = zeros(size(A_RM, 1), 1);
end


