function p = CIproj_p(c, delta, alphac, c_bd, deltastar_demean, Al, Au, sigma_l, sigma_u, deltaSigma, corr_m, corr_l, corr_u, eta, tol_r)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
lambda_l = Al*delta;
lambda_u = Au*delta;
lb = min(lambda_l);
ub = max(lambda_u);
mb = (lb + ub)/2;

deltastar = delta' + deltastar_demean;
lambdastar_l = deltastar*Al';
lambdastar_u = deltastar*Au';

% condition Delta
delta_sigma = sqrt(diag(deltaSigma));
maxdev = max(abs(deltastar_demean./delta_sigma'),[],2);
ind_Delta = (maxdev <= c_bd);

% condition proj
Tstar_l  = max(min((lambdastar_l - lb)./sigma_l',[],2), min((lb - lambdastar_u)./sigma_u',[],2));
Tstar_m  = max(min((lambdastar_l - mb)./sigma_l',[],2), min((mb - lambdastar_u)./sigma_u',[],2));
Tstar_u  = max(min((lambdastar_l - ub)./sigma_l',[],2), min((ub - lambdastar_u)./sigma_u',[],2));
ind_proj_l = (Tstar_l <= c); % not rejected
ind_proj_m = (Tstar_m <= c);
ind_proj_u = (Tstar_u <= c);

% condition c
% we only need to check condition c for select==1
%%%%%%% check this part
cLF = norminv(1 - eta);
selectl = (1 - ind_proj_l).*ind_Delta;
if sum(selectl) == 0
    ind_c_l = zeros(size(deltastar,1),1); 
else
    subl_deltastar = deltastar(selectl == 1,:);
    subl_lambdastar_l = lambdastar_l(selectl == 1,:);
    subl_lambdastar_u = lambdastar_u(selectl == 1,:);
    
    Tc_l1  = min((subl_lambdastar_l - lb)./sigma_l',[],2);
    Tc_l2  = min((lb - subl_lambdastar_u)./sigma_u',[],2);
    Tc_l   = max(Tc_l1, Tc_l2);
    [thl_1, thl_2] = CIcon_TNbounds(lb, subl_deltastar, Al, Au, sigma_l, sigma_u, corr_m, corr_l, corr_u, cLF, cLF, tol_r);
    phi_l = (normcdf(Tc_l) - normcdf(thl_1))./(normcdf(thl_2) - normcdf(thl_1));

    sub_ind_c_l = (phi_l <= 1 - alphac);
    ind_c_l = zeros(size(deltastar,1),1);
    ind_c_l(selectl == 1) = sub_ind_c_l;
end

selectm = (1 - ind_proj_m).*ind_Delta;
if sum(selectm) == 0
    ind_c_m = zeros(size(deltastar,1),1);
else
    subm_deltastar = deltastar(selectm == 1,:);
    subm_lambdastar_l = lambdastar_l(selectm == 1,:);
    subm_lambdastar_u = lambdastar_u(selectm == 1,:);
    
    Tc_m1  = min((subm_lambdastar_l - mb)./sigma_l',[],2);
    Tc_m2  = min((mb - subm_lambdastar_u)./sigma_u',[],2);
    Tc_m   = max(Tc_m1, Tc_m2);
    [thm_1, thm_2] = CIcon_TNbounds(mb, subm_deltastar, Al, Au, sigma_l, sigma_u, corr_m, corr_l, corr_u, cLF, cLF, tol_r);
    phi_m = (normcdf(Tc_m) - normcdf(thm_1))./(normcdf(thm_2) - normcdf(thm_1));

    sub_ind_c_m = (phi_m <= 1 - alphac);
    ind_c_m = zeros(size(deltastar,1),1);
    ind_c_m(selectm == 1) = sub_ind_c_m;
end

selectu = (1 - ind_proj_u).*ind_Delta;
if sum(selectu) == 0
    ind_c_u = zeros(size(deltastar,1),1);
else
    subu_deltastar = deltastar(selectu == 1,:);
    subu_lambdastar_l = lambdastar_l(selectu == 1,:);
    subu_lambdastar_u = lambdastar_u(selectu == 1,:);
    
    Tc_u1  = min((subu_lambdastar_l - ub)./sigma_l',[],2);
    Tc_u2  = min((ub - subu_lambdastar_u)./sigma_u',[],2);
    Tc_u   = max(Tc_u1, Tc_u2);
    [thu_1, thu_2] = CIcon_TNbounds(ub, subu_deltastar, Al, Au, sigma_l, sigma_u, corr_m, corr_l, corr_u, cLF, cLF, tol_r);
    phi_u = (normcdf(Tc_u) - normcdf(thu_1))./(normcdf(thu_2) - normcdf(thu_1));

    sub_ind_c_u = (phi_u <= 1 - alphac);
    ind_c_u = zeros(size(deltastar,1),1);
    ind_c_u(selectu == 1) = sub_ind_c_u;
end

ind_l = 1 - (1 - ind_c_l).*(1 - ind_proj_l);
ind_m = 1 - (1 - ind_c_m).*(1 - ind_proj_m);
ind_u = 1 - (1 - ind_c_u).*(1 - ind_proj_u);

p1  =  mean((1 - ind_l.*ind_m).*ind_Delta);
p2  =  mean((1 - ind_u.*ind_m).*ind_Delta);
p   =  max(p1, p2) + eta;
end

