function [sigma, dsigma] = fun_sigma(P)

brug = 1.5;
sigma = [P.sigma_neg*ones(P.n_neg+1,1)*(1-P.epsl_neg)^brug; zeros(P.n_sep-1,1); ...
    P.sigma_pos*ones(P.n_pos+1,1)*(1-P.epsl_pos)^brug];
dsigma = zeros(P.nj,1);