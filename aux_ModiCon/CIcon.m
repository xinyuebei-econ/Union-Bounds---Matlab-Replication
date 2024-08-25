function CIcon = CIcon(deltahat, deltaSigma, Al, Au, c_LF, alphac, tol, tol_r)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
kk = size(Al,1);
lambdaSigma = [Al; Au]*deltaSigma*[Al; Au]';
lambda_sigma = sqrt(diag(lambdaSigma));
corr_all = diag(lambda_sigma)^(-1)*lambdaSigma*diag(lambda_sigma)^(-1);
corr_m   = corr_all(1:kk, kk+1:end);
corr_l   = corr_all(1:kk, 1:kk);
corr_u   = corr_all(kk+1:end, kk+1:end);
sigma_l  = lambda_sigma(1:kk);
sigma_u  = lambda_sigma(kk+1:end);

lb = Al*deltahat - sigma_l*c_LF;
ub = Au*deltahat + sigma_u*c_LF;
lb = min(lb);
ub = max(ub);

mid = (lb + ub)/2;
lb_pt = min(Al*deltahat);
ub_pt = max(Au*deltahat);
% CI lower bound
rej = 1; theta = lb;
while rej == 1 && theta <= min(mid, lb_pt)  
    Tcl  = min((Al*deltahat - theta)./sigma_l);
    Tcu  = min((theta - Au*deltahat)./sigma_u);
    Tc = max( Tcl, Tcu);
    
    [th_1, th_2] = CIcon_TNbounds(theta, deltahat', Al, Au, sigma_l, sigma_u, corr_m, corr_l, corr_u, c_LF, c_LF, tol_r);
    t = norminv((1 - alphac)*normcdf(th_2) + alphac*normcdf(th_1));
    rej = (Tc > t);
    theta = theta + tol;
end
CIcon(1,1) = theta;

% CI upper bound
rej = 1; theta = ub;
while rej == 1 && theta >=  max(mid, ub_pt)    
    Tcl  = min((Al*deltahat - theta)./sigma_l);
    Tcu  = min((theta - Au*deltahat)./sigma_u);
    Tc = max( Tcl, Tcu);
    
    [th_1, th_2] = CIcon_TNbounds(theta, deltahat', Al, Au, sigma_l, sigma_u, corr_m, corr_l, corr_u, c_LF, c_LF, tol_r);
    t = norminv((1 - alphac)*normcdf(th_2) + alphac*normcdf(th_1));
    rej = (Tc > t);
    theta = theta - tol;
end
CIcon(1,2) = theta;

end

