function [delta0, deltaSigma, T_pre, Al, Au, index] = calibrateDGP(DGP, Type, M)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here
if DGP == 1
    load('./data/BCdata.mat','betahat','prePeriodIndices','sigma')
elseif DGP == 2
    load('./data/LWdata_EventStudy_Male.mat','betahat','prePeriodIndices','sigma')
elseif DGP == 3
    load('./data/LWdata_EventStudy_Female.mat','betahat','prePeriodIndices','sigma')
elseif DGP == 4
    load('./data/Christensen23_em.mat','betahat','prePeriodIndices','sigma')
elseif DGP == 5
    load('./data/Christensen23_match.mat','betahat','prePeriodIndices','sigma')
elseif DGP == 6
    load('./data/Dustmann_lemp.mat','betahat','prePeriodIndices','sigma')
end    
    
T_pre = length(prePeriodIndices);
deltaSigma = sigma(1:T_pre+1,1:T_pre+1);

ell = zeros( T_pre + 1, 1);% ell'*tau is the parameter of interest
ell(T_pre + 1) = 1; 
ell_post = ell(T_pre + 1:end,1);
gamma = sum(abs(ell_post'));
aux1 = -eye(T_pre);
aux2 =  eye(T_pre+1); aux2(1,:) = []; aux2(:,end) = [];
L_pre = aux1 + aux2;
A  = [ell'; [L_pre zeros(T_pre,1)]];    
Al = [ones(T_pre*2,1) gamma*M*[eye(T_pre); -eye(T_pre)]];
Au = [ones(T_pre*2,1) gamma*M*[eye(T_pre); -eye(T_pre)]]; 
Al = Al*A;
Au = Au*A;  
index = T_pre+1;
lambdasigma = sqrt(diag(Al*deltaSigma*Al'));

if Type == 1
    % parallel trend
    delta0 = zeros(T_pre+1,1);
elseif Type == 2
    % one large violation
    delta0 = zeros(T_pre+1,1);
    delta0(1) = 10*max(lambdasigma);
elseif Type == 3
    % small violations
    delta0 = [betahat(1:T_pre);0];
elseif Type == 4
    % linear trend
    delta0 = linspace(T_pre*10*max(lambdasigma),0,T_pre+1)';    
end

end

