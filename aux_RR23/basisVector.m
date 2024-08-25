function v = basisVector(index, size)
    if nargin < 2
        size = 1;
    end

    v = zeros(size, 1);
    v(index) = 1;
end

