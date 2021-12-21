
%this function produces samples from a multivariate conditional kernel density
%estimated distribution, assuming gaussian kernels are used to estimate the
%density.

%the response vector and the covariate vector can both be multivariate 
%can have variable bandwidths

%The method used here is based on the fact that the conditional density is
%a weighted mixture of normals.  A standard method for sampling from normal
%mixture is then used (sample probability of selecting a kernel based on
%weight, then sample from the corresponding distribution)

function [randsamp, w, distinfo, resind, ptindsel]=sampleCKDE_multcov_singres(datapts, bws, nreq, covind, covval) 

%INPUTS:
%datapts = [n x m] matrix of the sample data where n = no. data points, m = total number of
%variables 
%bws = [n x m] matrix of bandwidths for each data point 
%nreq = scalar indicating the no. of sample points required 
%covind = vector indicating index of covariates in 2nd dimension of datapts
%matrix.
%covval = [b x m] vector indicating the value of the covariate vector 

%ASSUMPTIONS:
%KDE is based on gaussian kernels
%product method has been used in the KDE (i.e. independent kernels for each
%dimension, which are then mulitplied to produce the density for each
%point)

npts = size(datapts,1);  %no. data points in sample set
nvars = size(datapts,2);  %total no. of variables 
resind = setdiff([1:nvars], covind);  %vector of the response variable 

%1.first determine the weight (i.e p(xi)/p(x))
[~, pxind] = MKDE_varbw_estdistn(datapts(:,covind), bws(:,covind), covval);
w = pxind./repmat(sum(pxind), size(pxind,1),1);   %weights on each of the normal components 

%2. now generate samples that indicate which component will be sampled
%from:

[~, tind] = sort(w);
distinfo = NaN*ones(npts, 3, size(covval,1));

for k = 1:size(distinfo,3)
    distinfo(:,:,k) = [w(tind(:,k),k), datapts(tind(:,k),resind), bws(tind(:,k), resind)];    
end
wmat = [zeros(1, size(covval,1)); cumsum(squeeze(distinfo(1:end-1,1,:))); ones(1, size(covval,1))];

%3. Now select points from the required components:
randsamp = NaN*ones(nreq, size(covval,1));
unos = rand(nreq, size(covval,1));
%determine which points should be taken for mean and variance of normal:
ptindsel = randsamp;
for k = 1:nreq    
    ptindsel(k,:) = sum((repmat(unos(k,:), size(wmat,1), 1) > wmat));        
end

for k = 1:size(covval,1)
    randsamp(:,k) = normrnd(distinfo(ptindsel(:,k),2,k), distinfo(ptindsel(:,k),3,k));
end


