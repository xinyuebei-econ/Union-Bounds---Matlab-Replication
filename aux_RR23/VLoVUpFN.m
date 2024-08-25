function VLoVUpVec = VLoVUpFN(eta, Sigma, A, b, z)
    c = LeeCFN(eta, Sigma);

    objective = (b - A * z) ./ (A * c);

    ACNegativeIndex = find((A * c) < 0);
    ACPositiveIndex = find((A * c) > 0);

    if isempty(ACNegativeIndex)
        VLo = -Inf;
    else
        VLo = max(objective(ACNegativeIndex));
    end

    if isempty(ACPositiveIndex)
        VUp = Inf;
    else
        VUp = min(objective(ACPositiveIndex));
    end
    VLoVUpVec = [VLo, VUp];
end

