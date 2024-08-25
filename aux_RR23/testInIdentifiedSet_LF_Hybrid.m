function reject = testInIdentifiedSet_LF_Hybrid(y, sigma, A, d, alpha, hybrid_list)
    sigmaTilde = sqrt(diag(A * sigma * A'));
    Atilde = (diag(1 ./ sigmaTilde)) * A;
    dtilde = (diag(1 ./ sigmaTilde)) * d;
    
    normalizedMoments = Atilde * y - dtilde;
    [maxMoment, maxLocation] = max(normalizedMoments);
    
    if maxMoment > hybrid_list.lf_cv
        reject = 1;
        return;
    else
        T_B = selectionMat(maxLocation, size(Atilde, 1), 'rows');
        iota = ones(size(Atilde, 1), 1);
        
        gamma = (T_B * Atilde)';
        Abar = Atilde - iota * (T_B * Atilde);
        dbar = (eye(size(dtilde, 1)) - iota * T_B) * dtilde;
        
        sigmabar = sqrt(gamma' * sigma * gamma);
        c = (sigma * gamma) / (gamma' * sigma * gamma);
        z = (eye(size(y, 1)) - c * gamma') * y;
        VLoVUpVec = VLoVUpFN(gamma, sigma, Abar, dbar, z);

        alphatilde = (alpha - hybrid_list.hybrid_kappa) / (1 - hybrid_list.hybrid_kappa);
        criticalVal = max(0, norminvp_generalized(1 - alphatilde, VLoVUpVec(1), VLoVUpVec(2), T_B * dtilde, sigmabar));
        reject = (maxMoment + T_B * dtilde > criticalVal);
    end
end

