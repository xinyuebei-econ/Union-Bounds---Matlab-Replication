function q = norminvp_generalized(p, l, u, mu, sd)
    if nargin < 5
        mu = 0;
        sd = 1;
    end

    lnormalized = (l - mu) / sd;
    unormalized = (u - mu) / sd;
    qnormalized = norminv((1-p)*normcdf(lnormalized) + p*normcdf(unormalized));
    q = mu + qnormalized * sd;
end

