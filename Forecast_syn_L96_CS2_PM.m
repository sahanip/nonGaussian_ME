%lead time forecasts Assimilation CS2 Proposed Method 

clear all
close all

%--------------------------------------------------------------------------
%first load and extract all the relevant information: 
load evaldata_syn_L96_CS2.mat %evaluation period data!
tstep = ts;  %can't use h because that's often used as a figure handle.

simlength = 100;  %first 1:simlength part of each IC is run!

x_obs_full = x_obs(:, 1:simlength+1,:);
x_true_full = x_true(:, 1:simlength+1,:);
errx_true_full = errx_true(:, 1:simlength+1,:);

clearvars -except tstep K Fcons hx hy obsfreq *_full noICs obserrvar simlength  

%now load the density estimation bandwidth data:
load bwresults_CS2_PM.mat

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX 
plotfigs = 0;  % plot if 1, don't plot if 0
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

maxleadtime = 15; 
nobins = 10;  %for the rankhistogram 

for ml = 1:noICs
    
    clearvars -except ml simlength maxleadtime *_ar bw bwm covindex datatrain plotfigs tstep K Fcons hx hy obsfreq *_full LSall CRPSall RMSEall ENSVARall ENSMEANERRall noICs obserrvar BINPOSall nobins *_anal
    
    %define the truth here:
    x_obs = x_obs_full(:,:,ml);
    x_true = x_true_full(:,:,ml);
    errx_true = errx_true_full(:,:,ml);

    %--------------------------------------------------------------------------
    %Input parameters:
    n = 1000;  % no. particles or ensemble members 
    x0var = 0; %error variance on x(0)
    novars = K;  %total no. of vars including the obs variable 
    obsind = [1 2 5 6];  %only every second variable observed
    hidind = [3 4 7 8 9];  %hidden states
    H = zeros(length(obsind), novars);  %define observation operator 
    for l = 1:length(obsind)
        H(l, obsind(l)) = 1;
    end
    dtrng = 0.05;
    rng2 = 0.5;
    varorder= [2:8];
    obserr = obserrvar(obsind,obsind);  %obs error statistics (mean and var) known perfectly and assume independence
    upstind = [1:K]; %index of internal states to update 
    zeroval = log(10^-8);  %this is for log score, when the estimated probability is zero.
    %--------------------------------------------------------------------------
    %Enkf parameters:
    locon = 0;
    tol = 10^-10;
    %--------------------------------------------------------------------------
    %STEP 1. First generate an ensemble of initial conditions, and preallocate
    %necessary matrices:
    x0_true = x_true(:,1);
    x0cov = x0var*eye(novars, novars);
    initstates = mvnrnd(repmat(x0_true', n,1), repmat(x0cov, [1, 1, n]))';

    %pre-allocate matrices:
    tempensstates = NaN*ones(novars, n, simlength+1);
    tempensstates(upstind,:,1) = initstates;
    ensstates = tempensstates; statesA = tempensstates; statesAact = tempensstates; statesAnew = statesA;
    Rall = NaN*ones(length(obsind), length(obsind), simlength);

    wghtall = NaN*ones(n, simlength+1);  %to store the weights 
    wghtall(:,1) = 1/n;  %all equally weighted first 

    %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    %state variable bounds - FOLLOWING WILL CHANGE DEPENDING ON THE PROBLEM
    xlb = -20; xub = 20;
    stbound =  xlb*ones(K, n);
    stbound(:,:,2) = xub*ones(K, n);
    stboundeq = NaN*ones(size(stbound));
    stdyncons = [];stdynconseq = [];
    xtbound = stbound(hidind,:,:);  %for the hidden states 
    ybound = stbound(obsind,:,:); %for the observed states
    %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    %--------------------------------------------------------------------------
    %Specify matrices for model error sampling:

    covind = [2];  %location of covariates in the data matrix input to CKDE 
    partchoice = 1;
    yevallength = 1000;
    yeval = linspace(-0.4,0.4,yevallength)';  %points at which to evaluate f(ey)
    f_ey = NaN*ones(yevallength, 4, simlength);  %to store f(ey) - ONLY GOING TO DO THIS FOR THE INDEP VARS
    eysample = NaN*ones(novars, n, simlength);  %to store the sampled error distribution 
    prodprobthresh = 10^-8;
    arvec = NaN*ones(novars, n,simlength);
    arvec(:,:,1) = repmat(errx_true(:,2), 1, n);
    errbnds = [min(datatrain(:,1))-0.05, max(datatrain(:,1))+0.05];  

    %--------------------------------------------------------------------------
    %STEP 2. Start the cycling process
    
    ls = NaN*ones(novars, simlength, maxleadtime);
    CRPS_calc = ls; 
    ensvar = ls;
    ensmeanerr = ls;
    binpos = ls;

    for i = 1:simlength 
        [ml,i]

        %--------------------------------------------------------------------------
        %STEP 1. Undertake one step of the model run with the ensemble from
        %previous time step:

        statesAnew(:,:,i) = statesA(:,:,i);

        sttemp = NaN*ones(K, n, obsfreq+1); sttemp(:,:,1) = statesAnew(:,:,i);  %at the last observed time
        for j = 1:obsfreq
            sttemp(:,:,j+1) = lorenz96(sttemp(:,:,j),Fcons,K, tstep,n);
        end
        tempensstates(:,:,i+1) = sttemp(:,:,end);

        %now check if any modelled ones are outside physical bounds + buffer
        keeppts = find(sum((and(tempensstates(:,:,i+1) > stbound(:,:,1), tempensstates(:,:,i+1) < stbound(:,:,2)))) == novars);
        unph = setdiff([1:n], keeppts);
        if length(keeppts) < n
            %redefine initial conditions:          
            statesAnew(:,unph,i) = statesAnew(:, keeppts(ceil(rand([1,length(unph)])*length(keeppts))), i);
    
            %run again:
            sttemp = NaN*ones(K, n, obsfreq+1); sttemp(:,:,1) = statesAnew(:,:,i);  %at the last observed time
            for j = 1:obsfreq
                sttemp(:,:,j+1) = lorenz96(sttemp(:,:,j),Fcons,K, tstep,n);
            end
            tempensstates(:,:,i+1) = sttemp(:,:,end);
        end

        %--------------------------------------------------------------------------
        %STEP 2. Determine covariates for conditional distribution for each variable:

        fixedres = 1;
        
        for p = 1:novars
            if p == 1
                pp = novars;
            else
                pp = p - 1;
            end
            
            %now sample values:
            clear temperrest 
            
            arvec(p,:,i) = statesA(p,:,i) - tempensstates(p,:,i);
            arvec(p,:,i) = max([arvec(p,:,i); errbnds(1)*ones(1,n)]);
            arvec(p,:,i) = min([arvec(p,:,i); errbnds(2)*ones(1,n)]);
                        
            covvaltemp = [statesAnew(p,:,i)', statesAnew(pp,:,i)',arvec(p,:,i)'];
            temperrest = sampleCKDE_multcov_singres(datatrain(:,[fixedres,covindex(fixedres,1:size(covindex,2))]), repmat(bw(fixedres,:), size(datatrain,1),1), 1, covindex, covvaltemp); 
            
            if isempty(temperrest)                
                %redefine covariates: locate closest point in covariate space and then sample:
                [~,mc] = min(sum((datatrain(:, covindex) - repmat(covvaltemp, size(datatrain,1),1)).^2,2));
                covvaltemp = datatrain(mc,covindex);
                temperrest = sampleCKDE_multcov_singres(datatrain(:,[fixedres,covindex(fixedres,1:size(covindex,2))]), repmat(bw(fixedres,:), size(datatrain,1),1), 1, covindex, covvaltemp);             
            end

            %add error and check if physical:
            %Don't worry about physicality of error value just yet, leave this to longer lead times 
            eysample(p,:,i) = temperrest;
            ensstates(p,:,i+1) = eysample(p,:,i) + tempensstates(p,:,i+1); %range of q based on model error for q
            unph = find(or(ensstates(p,:,i+1) < stbound(p,:,1), ensstates(p,:,i+1) > stbound(p,:,2)));
            while ~isempty(unph)
                %resample:
                eysample(p,unph,i) = sampleCKDE_multcov_singres(datatrain(:,[fixedres,covindex(fixedres,1:size(covindex,2))]), repmat(bw(fixedres,:), size(datatrain,1),1), 1, covindex, covvaltemp(unph,:));             
                ensstates(p,unph,i+1) = eysample(p,unph,i) + tempensstates(p,unph,i+1); %range of q based on model error for q
                unph = find(or(ensstates(p,:,i+1) < stbound(p,:,1), ensstates(p,:,i+1) > stbound(p,:,2)));
            end 
        end
        
        %--------------------------------------------------------------------------
        %STEP 4. Undertake DA Update 

        %USE ETKF:
        Rall(:,:,i) = obserr;
        lambdanew = 1;  %no inflation 
        CLOC = ones(novars, novars); %no localization 
        [statesA(:,:,i+1), ~, ~] = ETKF_inf_loc(x_obs(obsind,i+1), ensstates(:,:,i+1), Rall(:,:,i), n,H, tol, lambdanew, CLOC);
        wghtall(:,i+1) = 1/n;  %equally weighted for enkf  

        %--------------------------------------------------------------------------
        %STEP 6. MODEL FORWARD IN TIME FOR MAX NO. OF LEAD TIMES &
        %CALCULATE SCORES:
        
        %first determine max lead time:
        maxleadtime2 = min(maxleadtime, simlength-i);
        
        for lt = 1:maxleadtime2 
            
            clear xfit ffit mtemp etemp predic2 artemp arcov prediccov 
            if lt == 1
                %already have one step ahead prediction:
                predic = ensstates(:,:,i+1);
                eotemp = (statesA(:,:,i+1) - tempensstates(:,:,i+1));                                
            else
                %generate long leadtime forecast:    
                
                %first determine covariate (since predic will be
                %overwritten): 
                prediccov = predic;
                arcov = eotemp;
                
                %now generate model simulation:
                sttemp = NaN*ones(K, n, obsfreq+1); sttemp(:,:,1) = predic;  %at the last observed time
                for j = 1:obsfreq
                    sttemp(:,:,j+1) = lorenz96(sttemp(:,:,j),Fcons,K, tstep,n);
                end
                predic = sttemp(:,:,end);
                
                %now sample errors:
                fixedres = 1;
                for p = 1:novars
                    if p == 1
                        pp = novars;
                    else
                        pp = p - 1;
                    end

                    %now sample values:
                    clear temperrest 

                    %now redefine error region if outside: 
                    arcov(p,:) = max([arcov(p,:); errbnds(1)*ones(1,n)]);
                    arcov(p,:) = min([arcov(p,:); errbnds(2)*ones(1,n)]);
                        
                    covvaltemp = [prediccov(p,:)', prediccov(pp,:)',arcov(p,:)'];
            
                    temperrest = sampleCKDE_multcov_singres(datatrain(:,[fixedres,covindex(fixedres,1:size(covindex,2))]), repmat(bw(fixedres,:), size(datatrain,1),1), 1, covindex, covvaltemp); 
            
                    if isempty(temperrest)
                        %redefine covariates: locate closest point in covariate space and then sample:
                        [~,mc] = min(sum((datatrain(:, covindex) - repmat(covvaltemp, size(datatrain,1),1)).^2,2));
                        covvaltemp = datatrain(mc,covindex);
                        temperrest = sampleCKDE_multcov_singres(datatrain(:,[fixedres,covindex(fixedres,1:size(covindex,2))]), repmat(bw(fixedres,:), size(datatrain,1),1), 1, covindex, covvaltemp);             
                    end

                    %add error and check if physical:
                    %Don't worry about physicality of error value just yet, leave this to longer lead times 
                    eotemp(p,:) = temperrest;
                    predic2(p,:) = eotemp(p,:) + predic(p,:); %range of q based on model error for q
                    unph = find(or(predic2(p,:) < stbound(p,:,1), predic2(p,:) > stbound(p,:,2)));
                    while ~isempty(unph)
                        %resample:
                        eotemp(p,unph) = sampleCKDE_multcov_singres(datatrain(:,[fixedres,covindex(fixedres,1:size(covindex,2))]), repmat(bw(fixedres,:), size(datatrain,1),1), 1, covindex, covvaltemp(unph,:));                             
                        predic2(p,unph) = eotemp(p,unph) + predic(p,unph); %range of q based on model error for q
                        unph = find(or(predic2(p,:) < stbound(p,:,1), predic2(p,:) > stbound(p,:,2)));
                    end 
                end

                %checks to make sure covering everything 
                if ~isempty(find(isnan(eotemp)))
                    error('NAN ERRORS FOUND!')
                end
                
                %NO UPDATING:
                predic = predic2; 
                clear predic2   
            end

            %now calculate scores:
                       
            %LOG SCORE (NEGATIVE VERSION) & CRPS:
            zeroval = log(10^-8);  %this is for log score, when the estimated probability is zero.
            for k = 1:novars
                clear xfit ffit 
                [ffit, xfit] = ksdensity(squeeze(predic(upstind(k),:)), 'Weights', wghtall(:,i));   
                [ls(k,i,lt), ~] = logscore(xfit', ffit', x_true(k,i+lt), zeroval); %Note, calculating on the truth...
                [CRPS_calc(k,i,lt)] = CRPS(predic(k,:)', wghtall(:,i), x_true(k,i+lt));
            end
            ensvar(:,i,lt) = var(predic,[],2);  %ensemble variance
            ensmeanerr(:,i,lt) = mean(predic,2) - x_true(:,i+lt);  %same error in ensemble mean to calculate RMSE later!  
            
            
        end
    end     
        
       
    %now collate all stats:  
    
    %MAKING IT NEGATIVE:
    LSall{ml} = -ls;
    
    %CRPS:
    CRPSall{ml} = CRPS_calc;
    
    %ENS MEAN::
    ENSMEANERRall{ml} = ensmeanerr;
    
    %ENS var:    
    ENSVARall{ml} = ensvar; 
    
    %analysis:
    ENSMEANERRall_anal{ml} = squeeze(mean(statesA,2)) - x_true;
    
    ENSVARall_anal{ml} = squeeze(var(statesA,[],2));
end

save('forecastresults_CS2_PM.mat', 'LSall', 'CRPSall', 'ENSMEANERRall', 'ENSVARall', 'maxleadtime', 'nobins', '*_anal', 'x_true')

