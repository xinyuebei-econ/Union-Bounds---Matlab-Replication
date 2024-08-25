function [critical_value] = compute_least_favorable_cv(X_T, sigma, hybrid_kappa, sims, rowsForARP)
    if nargin < 4
        sims = 1000;
    end
    
    if nargin < 5
        rowsForARP = [];
    end
    
    if ~isempty(rowsForARP)
        if isvector(X_T)
            X_T = X_T(rowsForARP);
        else
            X_T = X_T(rowsForARP, :);
        end
        sigma = sigma(rowsForARP, rowsForARP);
    end
    
    compute_eta = @(b, f, C) compute_eta_lp(b, f, C);
    
    rng(0);
    
    if isempty(X_T)
        xi_draws = mvnrnd(zeros(sims, size(sigma, 1)), sigma) ./ sqrt(diag(sigma))';
        eta_vec = max(xi_draws, [], 1);
        critical_value = quantile(eta_vec, 1 - hybrid_kappa);
    else
        sdVec = sqrt(diag(sigma));
        dimDelta = size(X_T, 2);
        f = [1; zeros(dimDelta, 1)];
        C = -[sdVec, X_T];
        xi_draws = mvnrnd(zeros(sims, size(sigma, 1)), sigma)';
        
        eta_vec = zeros(sims, 1);
        for i = 1:sims
            eta_vec(i) = compute_eta(-xi_draws(:, i), f, C);
        end
        
        critical_value = quantile(eta_vec, 1 - hybrid_kappa);
    end
end

function eta = compute_eta_lp(b, f, C)
%     The problem is min f*x such that C*x<=b
%     linprog = lpSolveAPI.make.lp(0, length(f));
%     lpSolveAPI.set.objfn(linprog, f); f is the objective function 
%     for r = 1:size(C, 1)
%         lpSolveAPI.add.constraint(linprog, C(r, :), '<=', b(r));
%     end

%     lpSolveAPI.set.bounds(linprog, -Inf(length(f), 1), 1:length(f));
%     lpSolveAPI.lp.control(linprog, 'min', 'dual', 'dantzig', 'neutral', 10);
%     error_flag = lpSolveAPI.solve.lpExtPtr(linprog);
%     if error_flag
%         eta = NaN;
%     else
%         eta = lpSolveAPI.get.objective(linprog);
%     end
    options = optimoptions('linprog','Display','off');
    save('lp_compute_least_favorable_cv','f','C','b')
    stop
    [~, eta] = linprog(f, C, b,[],[],[],[], options);
end

