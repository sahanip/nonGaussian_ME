
function fcond_est = CKDE_varbw_estdistn(X, bw, bwm, xp, covind)

%ASSUMPTIONS:
%Normal kernel 
%normalise X for determining k nearest neighbour 
%this one only looks at the covariate dimensions to determine the nearest
%neighbours.

%INPUTS:
%X = n x m matrix of sample data where n = no. of data points, m = no. of variables 
%bw = n x m matrix of bandwidths for each sample point for the joint distn in numerator
%bwm = n x l matrix of bandwidths for each sample point for the marginal distn in denominator 
%xp = yn x m matrix of points at which to estimate density 
%nostds = no. of standard deviations to equate the distance to the farthest
%k nearest neighbour.
%covind = vector indicating columns that are covariates.  Must be such that covind corresponds to the columns in bwm 

n = size(X,1);  %no. of sample points 
v = size(X,2);  %no. of dimensions 
m = length(covind);  %no. of covariates
yn = size(xp,1);   %no. of points to evaluate distn 

indepkern = NaN*ones(n, v, yn);  
prodkern = NaN*ones(n,yn);
f_est = NaN*ones(yn,1);

%for the numerator:
for j = 1:yn
    indepkern(:,:,j) = normpdf(repmat(xp(j,:), n, 1),X,bw);
    prodkern(:,j) = prod(indepkern(:,:,j),2);
    f_est(j,1) = mean(prodkern(:,j));  
end

%for the denominator:
indepkernden = NaN*ones(n, m, yn);  
prodkernden = NaN*ones(n,yn);
f_estden = NaN*ones(yn,1);

for j = 1:yn
    indepkernden(:,:,j) = normpdf(repmat(xp(j,covind), n, 1),X(:,covind),bwm);
    prodkernden(:,j) = prod(indepkernden(:,:,j),2);
    f_estden(j,1) = mean(prodkernden(:,j));  
end

fcond_est = f_est./f_estden;   








    




