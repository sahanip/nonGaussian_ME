function[xnew, enew] = resample_region_v2(i, p, xtest, eysample, datatrain, bwadapt,covindex, fixedres, ballsize,pdensthresh,modval) 

%this script resamples parameters whenever they are outside a desired
%region.  Points are mapped to nearest training data point when this
%occurs. 

%same as resample_region.m except that it assumes the covariates
%are only x_k and e_x_k (i.e. 2 covs only) and nsamp points is reduced to
%500

%INPUTS:
%i = time index
%p = variable of interest 
%pp = Xk-1 
%xtest = k x t vector of states 
%eysample = k x t vector of errors 
%ballsize = max distance in each direction to search 
%modval = sttemp(p,1,end), for example.

nsamp = 500;

disttemp = abs(datatrain(:,covindex)  - repmat([xtest(p,i+1), eysample(p,i+1)], size(datatrain,1),1));
pdens = length(find(((disttemp(:,1) < ballsize(1)) + (disttemp(:,2) < ballsize(2))) == 2));
if pdens < pdensthresh
    %sample more points and check:
    tmppts = sampleCKDE(datatrain(:,[fixedres,covindex]), bwadapt, nsamp, covindex, [xtest(p,i), eysample(p,i)], 0,0);
    xtmppts = modval + tmppts;
    clear pdens
    for lm = 1:nsamp
        disttemp = abs(datatrain(:,covindex)  - repmat([xtmppts(lm), tmppts(lm)], size(datatrain,1),1));
        pdens(lm) = length(find(((disttemp(:,1) < ballsize(1)) + (disttemp(:,2) < ballsize(2))) == 2));
    end
    pcheck = find(pdens > pdensthresh);
    if isempty(pcheck)
        %now map onto training data:
        [~, closestpt] = min(sum((datatrain(:,covindex)  - repmat([xtest(p,i), eysample(p,i)], size(datatrain,1),1)).^2,2));
        xtest(p,i+1) = datatrain(closestpt(1),covindex(1));
        eysample(p,i+1) = datatrain(closestpt(1), covindex(2));
    else
        ps = randsample([1:length(pcheck)],1);
        xtest(p,i+1) = xtmppts(pcheck(ps));
        eysample(p,i+1) = tmppts(pcheck(ps));
    end
end

xnew = xtest(p,i+1);
enew = eysample(p,i+1);
