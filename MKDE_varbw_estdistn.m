
function [f_est,prodkern] = MKDE_varbw_estdistn(X, bw, xp)

%ASSUMPTIONS:
%Normal kernel 
%normalise X for determining k nearest neighbour 
%this one only looks at the covariate dimensions to determine the nearest
%neighbours.

%INPUTS:
%X = n x m matrix of sample data where n = no. of data points, m = no. of variables 
%bw = n x m matrix of bandwidths for each sample point 
%xp = yn x m matrix of points at which to estimate density 
%nostds = no. of standard deviations to equate the distance to the farthest
%k nearest neighbour.
%covinds = vector indicating columns that are covariates 

n = size(X,1);  %no. of sample points 
v = size(X,2);  %no. of dimensions 
yn = size(xp,1);   %no. of points to evaluate distn 

%If using a normal kernel, bandwidth for dimension x1 is the standard
%deviation of the independent univariate normal kernel

indepkern = NaN*ones(n, v, yn);  
prodkern = NaN*ones(n,yn);
f_est = NaN*ones(yn,1);

for j = 1:yn
    indepkern(:,:,j) = normpdf(repmat(xp(j,:), n, 1),X,bw);
    prodkern(:,j) = prod(indepkern(:,:,j),2);
    f_est(j,1) = mean(prodkern(:,j));  
end











    




