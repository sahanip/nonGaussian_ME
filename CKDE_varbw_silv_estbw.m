
%This estimates bandwidths scaled by local factors based on the method in
%page 101 of Silverman's Density Estimation for Statistical and Data
%analysis

function [bw,bwm, bwfactors] = CKDE_varbw_silv_estbw(X, covind, initbw, initbwm, alpha)

%ASSUMPTIONS:
%Normal kernel 
%h is not re-optimised.
%initial bandwidth is fixed (not adaptive)

%INPUTS:
%X = n x m matrix of sample data where n = no. of data points, m = no. of variables 
%initbw = 1 x m vector of initial bandwidth estimate (fixed bandwidth) 
%initbwm = 1 x v vector of initial bandwidth estimate for covariates only  (fixed bandwidth) 
%alpha = scalar sensitivity parameter, usually 1/2
%vector indicating position of covariates.

n = size(X,1);  %no. of sample points 
v = size(X,2);  %no. of dimensions (total)
v2 = length(covind);  %no. covariates 

% %--------------------------------------------------------------------------
%first calculate the value of g:
%need to get density estimates at points of interest:

bwmat = repmat(initbw, n, 1);
bwm_mat = repmat(initbwm, n, 1);
f_estinit = CKDE_varbw_estdistn(X, bwmat, bwm_mat, X, covind);

g = exp(mean(log(f_estinit)));
lambda = (f_estinit*(1/g)).^-alpha;

%bwfactors = repmat(lambda, 1, v);
bwfactors = lambda;
bw = repmat(bwfactors, 1, v).*repmat(initbw, n,1);
bwm = repmat(bwfactors, 1, v2).*repmat(initbwm, n,1);

