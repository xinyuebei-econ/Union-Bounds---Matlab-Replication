function rounded = roundeps(x, eps)
    if nargin < 2
        eps = 10^(-10);
    end

    if abs(x) < eps
        rounded = 0;
    else
        rounded = x;
    end
end
