function[xnew, enew] = resample_region(i, p, pp, xtest, eysample, datatrain, bwadapt,covindex, fixedres, ballsize,pdensthresh,modval) 

%this script resamples parameters whenever they are outside a desired
%region.  Points are mapped to nearest training data point when this
%occurs. 

%ASSUMES A SPECIFIC COVARIATE STRUCTURE - XK, XK-1, EXK

%INPUTS:
%i = time index
%p = variable of interest 
%pp = Xk-1 
%xtest = k x t vector of states 
%eysample = k x t vector of errors 
%ballsize = max distance in each direction to search 
%modval = sttemp(p,1,end), for example.

%to decide on the Xk-1 
if p > 1  %following is for the xtest(pp) covariate, since this doesn't exist yet for p = 1, use value at previous time as a "good enough" estimate. 
    j = i+1;
else
    j = i;
end

nsamp = 1000;

disttemp = abs(datatrain(:,covindex)  - repmat([xtest(p,i+1), xtest(pp,j), eysample(p,i+1)], size(datatrain,1),1));
pdens = length(find(((disttemp(:,1) < ballsize(1)) + (disttemp(:,2) < ballsize(2)) + (disttemp(:,3) < ballsize(3))) == 3));
if pdens < pdensthresh
    %sample more points and check:
    tmppts = sampleCKDE(datatrain(:,[fixedres,covindex]), bwadapt, nsamp, covindex, [xtest(p,i), xtest(pp,i), eysample(p,i)], 0,0);
    xtmppts = modval + tmppts;
    clear pdens
    for lm = 1:nsamp
        disttemp = abs(datatrain(:,covindex)  - repmat([xtmppts(lm), xtest(pp,j), tmppts(lm)], size(datatrain,1),1));
        pdens(lm) = length(find(((disttemp(:,1) < ballsize(1)) + (disttemp(:,2) < ballsize(2)) + (disttemp(:,3) < ballsize(3))) == 3));
    end
    pcheck = find(pdens > pdensthresh);
    if isempty(pcheck)
        %now map onto training data:
        [~, closestpt] = min(sum((datatrain(:,covindex(1:2))  - repmat([modval, xtest(pp,j)], size(datatrain,1),1)).^2,2));
        xtest(p,i+1) = datatrain(closestpt(1),covindex(1));
        eysample(p,i+1) = datatrain(closestpt(1), covindex(3));
    else
        ps = randsample([1:length(pcheck)],1);
        xtest(p,i+1) = xtmppts(pcheck(ps));
        eysample(p,i+1) = tmppts(pcheck(ps));
    end
end

xnew = xtest(p,i+1);
enew = eysample(p,i+1);
