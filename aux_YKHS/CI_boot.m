function [CI_YKHS, CI_Boot] = CI_boot(deltahat, data, Sigma, A_l, A_u, alpha, n, m, B)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here  
lambdahatl = A_l*deltahat;
lambdahatu = A_u*deltahat;
dmin = (1-sqrt(m)/sqrt(n))*(min(lambdahatl) - lambdahatl);
dmax = (1-sqrt(m)/sqrt(n))*(max(lambdahatu) - lambdahatu);

rng('default')
select = datasample(1:n,m,'Replace',false);
delta_hat_m =  mean(data(select,:))';

lambda_l_m = min(A_l*delta_hat_m);
lambda_u_m = max(A_u*delta_hat_m);

delta_boot  = mvnrnd(deltahat, Sigma, B);
lambda_bootl = min(delta_boot*A_l' + dmin',[],2);
lambda_bootu = max(delta_boot*A_u' + dmax',[],2);
    
CI_YKHS = [lambda_l_m - sqrt(n/m)*quantile(lambda_bootl - min(lambdahatl), 1 - alpha/2), ...
           lambda_u_m - sqrt(n/m)*quantile(lambda_bootu - max(lambdahatu), alpha/2)];
CI_Boot = [min(lambdahatl) - quantile((lambda_bootl - min(lambdahatl)), 1 - alpha/2), ...
           max(lambdahatu) - quantile((lambda_bootu - max(lambdahatu)), alpha/2)];
end

