function [th_1, th_2] = CIcon_TNbounds(theta, deltahat, Al, Au, sigma_l, sigma_u, corr_m, corr_l, corr_u, cLFl, cLFu, tol_r)
%UNTITLED3 Summary of this function goes here
%   sigma_l,sigma_u are std of lambdahat_l and lambdahat_u
%%
B = size(deltahat,1);
lambdahat_l = deltahat*Al';
lambdahat_u = deltahat*Au';

Tlb = (lambdahat_l - theta)./sigma_l';
Tub = (theta - lambdahat_u)./sigma_u';
[Tl, bl] = min(Tlb, [], 2);
[Tu, bu] = min(Tub, [], 2);

corr_m_blb = corr_m(bl,:);
corr_m_bub = corr_m(:,bu)';
%%
% for Tl > Tu 
% TS
tTS_l   = 10^10*ones(B,size(Au,1));
tTS_ll  = (1 + corr_m_blb).^(-1).*(Tub + corr_m_blb.*Tl);
tTS_l(1 + corr_m_blb > tol_r) = tTS_ll(1 + corr_m_blb > tol_r);
    
th_1_1 = min(tTS_l,[],2).*(Tl >= Tu);
th_1_1(th_1_1 == 10^10) = -10^10;
th_1_1 = th_1_1.*(Tl >= Tu);

% LF
tLFl = cLFl.*ones(B,1);

% B
taux_1 = (1 - corr_l(bl,:)).^(-1).*(Tlb - corr_l(bl,:).*Tl);
tB_l2  = 10^10*ones(B,size(Al,1));
tB_l2(1 >  corr_l(bl,:) + tol_r)  = taux_1(1 >  corr_l(bl,:) + tol_r); 

th_2_1 = min([tLFl, tB_l2],[],2).*(Tl >=  Tu);
%%
% for the upper bound
% TS
tTS_u   =  10^10*ones(B,size(Al,1));
tTS_u1  = (1 + corr_m_bub).^(-1).*(Tlb + corr_m_bub.*Tu);
tTS_u(1 + corr_m_bub >  tol_r) = tTS_u1(1 + corr_m_bub >  tol_r);

th_1_2 = min(tTS_u,[],2).*(Tl < Tu);
th_1_2(th_1_2 == 10^10) = -10^10;
th_1_2 = th_1_2.*(Tl < Tu);

% LF
tLFu = cLFu.*ones(B,1);

% B
taux_2 = (1 - corr_u(bu,:)).^(-1).*(Tub - corr_u(bu,:).*Tu);
tB_u2  = 10^10*ones(B,size(Au,1));
tB_u2(1 >  corr_u(bu,:) + tol_r) = taux_2(1 >  corr_u(bu,:) + tol_r);

th_2_2 = min([tLFu, tB_u2],[],2).*(Tl < Tu);

th_1 = th_1_1 + th_1_2;
th_2 = th_2_1 + th_2_2;
end

