function resultsGrid = testOverThetaGrid(betahat, sigma, A, d, thetaGrid, numPrePeriods, alpha, testFn, hybrid_list)
    % Tests whether values in a grid lie in the identified set.

    testTauInSetFn = @(theta) testFn(betahat - basisVector(numPrePeriods + 1, length(betahat)) * theta, sigma, A, d, alpha, hybrid_list);
    testResultsGrid = nan(1,length(thetaGrid)); % results is reject
    for i = 1:length(thetaGrid)
        testResultsGrid(1,i) = testTauInSetFn(thetaGrid(i));
    end
    testValsGrid = thetaGrid;
    resultsGrid = [testValsGrid', testResultsGrid'];
end


