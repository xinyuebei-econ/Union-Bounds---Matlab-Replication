function [CI_set, CI_pt] = CI_boot_pt(deltahat, data, Sigma, A_l, A_u, alpha, n, m, B)
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

delta_boot   = mvnrnd(deltahat, Sigma, B);
lambda_bootl = min(delta_boot*A_l' + dmin',[],2);
lambda_bootu = max(delta_boot*A_u' + dmax',[],2);
    
CI_set = [lambda_l_m - sqrt(n/m)*quantile(lambda_bootl - min(lambdahatl), 1 - alpha/2), ...
           lambda_u_m - sqrt(n/m)*quantile(lambda_bootu - max(lambdahatu), alpha/2)];
       
omegahat =   lambda_u_m - sqrt(n/m)*quantile(lambda_bootu - max(lambdahatu),1/2) ...
           - lambda_l_m + sqrt(n/m)*quantile(lambda_bootl - min(lambdahatl),1/2);
omegahat = max(omegahat,0);
rho = (log(m))^(-1)*sqrt(m/n)/max(quantile(lambda_bootu,3/4) - quantile(lambda_bootu,1/4), ...
                                  quantile(lambda_bootl,3/4) - quantile(lambda_bootl,1/4)); 
p = 1 - normcdf(rho*omegahat)*alpha;

CI_pt = [lambda_l_m - sqrt(n/m)*quantile(lambda_bootl - min(lambdahatl), p), ...
           lambda_u_m - sqrt(n/m)*quantile(lambda_bootu - max(lambdahatu), 1 - p)];
end

