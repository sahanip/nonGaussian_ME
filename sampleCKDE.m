%this function produces samples from a multivariate conditional kernel density
%estimated distribution, assuming gaussian kernels are used to estimate the
%density.

%the response vector and the covariate vector can both be multivariate 
%can have variable bandwidths

%The method used here is based on the fact that the conditional density is
%a weighted mixture of normals.  A standard method for sampling from normal
%mixture is then used (sample probability of selecting a kernel based on
%weight, then sample from the corresponding distribution)

function [randsamp, w, distinfo, resind, ftest]=sampleCKDE(datapts, bws, nreq, covind, covval, chk, xtest) 

%INPUTS:
%datapts = [n x m] matrix of the sample data where n = no. data points, m = total number of
%variables 
%bws = [n x m] matrix of bandwidths for each data point 
%nreq = scalar indicating the no. of sample points required 
%covind = vector indicating index of covariates in 2nd dimension of datapts
%matrix.
%covval = [1 x m] vector indicating the value of the covariate vector 
%chk = 1 if you want to plot the mixture distribution from which we are
%sampling, zero else.  Can only do this if the response variable is
%univariate. 
%xtest must be specified if chk == 1, for plotting of the distribution 

%ASSUMPTIONS:
%KDE is based on gaussian kernels
%product method has been used in the KDE (i.e. independent kernels for each
%dimension, which are then mulitplied to produce the density for each
%point)

nvars = size(datapts,2);  %total no. of variables 
resind = setdiff([1:nvars], covind);  %vector of the response variable 

%1.first determine the weight (i.e p(xi)/p(x))
[~, pxind] = MKDE_varbw_estdistn(datapts(:,covind), bws(:,covind), covval);
w = pxind/sum(pxind);   %weights on each of the normal components 

%2. now generate samples that indicate which component will be sampled
%from:
unos = rand(nreq,1);
[~, tind] = sort(w);
for k = 1:length(resind)
    distinfo(:,:,k) = [w(tind), datapts(tind,resind(k)), bws(tind, resind(k))];
end
wmat = [0; cumsum(distinfo(1:end-1,1,1)); 1];

%3. Now select points from the required components:
randsamp = [];
for k = 1:size(distinfo,1)
    sampno = length(find(and(unos > wmat(k), unos <= wmat(k+1))));
    if sampno > 0 
        
        randsamp = [randsamp; mvnrnd(repmat(squeeze(distinfo(k,2,:))', sampno, 1),(diag(squeeze(distinfo(k,3,:))')).^2)];   
    end
end

if and(resind == 1, chk == 1)
    %if response variable is univariate
    %calculate the pdf over a reasonable set of points:
    mutest = squeeze(distinfo(:,2,:));
    sigmatest = reshape(squeeze(distinfo(:,3,:)), [1,1,size(distinfo,1)]);
    
    for k = 1:length(xtest)        
        ftemp = mvnpdf(xtest(k)*ones(size(distinfo,1),1), mutest, sigmatest.^2);  
        ftest(k) = sum(distinfo(:,1).*ftemp);
    end
else
    ftest = [];  
end

    
    



