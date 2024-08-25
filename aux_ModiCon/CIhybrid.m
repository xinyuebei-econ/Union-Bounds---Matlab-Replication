function [CI_h, CI_c, CI_p, delta1, c] = CIhybrid(deltahat, deltaSigma, Al, Au, alpha, alphac, eta, B, Blarge, tol, tol_r, index)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
rng(0,'twister')
gs = GlobalSearch('Display','off','NumTrialPoints',3000);
deltastar_demean_large = mvnrnd(zeros(size(deltahat)), deltaSigma, Blarge);
deviation = max(abs(deltastar_demean_large./sqrt(diag(deltaSigma))'),[],2);
c_bd = quantile(deviation, 1 - eta);
sigma_delta = sqrt(diag(deltaSigma));        
lambdaSigma = [Al; Au]*deltaSigma*[Al; Au]';
lambdasigma = sqrt(diag(lambdaSigma));
kk = size(Al,1);
sigma_l  = lambdasigma(1:kk);
sigma_u  = lambdasigma(kk+1:end);
corr_all = diag(sqrt(diag(lambdaSigma)))^(-1)*lambdaSigma*diag(sqrt(diag(lambdaSigma)))^(-1);
corr_m   = corr_all(1:kk, kk+1:end);
corr_l   = corr_all(1:kk, 1:kk);
corr_u   = corr_all(kk+1:end, kk+1:end);

lb  = deltahat - sigma_delta*c_bd;
ub  = deltahat + sigma_delta*c_bd;
delta1 = deltahat;

lb(index) = 0;            ub(index) = 0;            delta1(index) = 0;

cl = 0; cu = norminv(1 - alpha/2); c  = (cl + cu)/2;
delta_fea = [];
c_fea     = [];
obj_large = @(delta, c_check) (alpha - CIproj_p(c_check, delta, alphac, c_bd,...
    deltastar_demean_large, Al, Au, sigma_l, sigma_u, deltaSigma, corr_m, corr_l, corr_u, eta, tol_r))*100;

k = 1;
while cu - cl > tol  
    rng(k,'twister')
    k = k+1;

    deltastar_demean = mvnrnd(zeros(size(deltahat)), deltaSigma, B);
    obj = @(delta) (alpha - CIproj_p(c, delta, alphac, c_bd, deltastar_demean, Al, Au, sigma_l,...
        sigma_u, deltaSigma, corr_m, corr_l, corr_u, eta, tol_r))*100;   
    p1  = obj_large(delta1, c);
    if  p1 >= 0
        problem = createOptimProblem('fmincon','x0', delta1,'objective', obj, 'lb', lb, 'ub', ub);  
        [ delta1, p] = run(gs, problem);
        p1 = obj_large(delta1,c);
        if p1 >= 0 
            cu = c;
        else 
            cl = c;
            delta_fea = [delta_fea delta1];
            c_fea     = [c c_fea];
        end
    else
        cl = c;
        delta_fea = [delta_fea delta1];
        c_fea     = [c c_fea];
    end   
    c = (cl + cu)/2;    
end

cl = c; cu = norminv(1-alpha/2);
while cu - cl > tol
    c = (cl + cu)/2;
    p = obj_large(delta1, c);
    if p >= 0
        cu = c;
    elseif p < 0
        cl = c;
    end
end            
lambdahat_l = Al*deltahat; lambdahat_u = Au*deltahat;

c_LF  = norminv(1 - eta);
CI_p  = [min(lambdahat_l - c*sigma_l) max(lambdahat_u + c*sigma_u)];
CI_c  = CIcon(deltahat, deltaSigma, Al, Au, c_LF, alphac, tol, tol_r);
CI_h  = [min(CI_p(1),CI_c(1)) max(CI_p(2),CI_c(2))];
end

