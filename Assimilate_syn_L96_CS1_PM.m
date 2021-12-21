%CS1 proposed method assimilation 

clear all
close all

%--------------------------------------------------------------------------
%first load and extract all the relevant information:
load evaldata_syn_L96_CS1.mat  %evaluation period data!
tstep = ts;  %can't use h because that's often used as a figure handle.

simlength = 100;  %first 1:simlength part of each IC is run!

x_obs_full = x_obs(:, 1:simlength+1,:);
x_true_full = x_true(:, 1:simlength+1,:);
errx_true_full = errx_true(:, 1:simlength+1,:);

clearvars -except tstep K Fcons hx hy obsfreq *_full noICs obserrvar simlength   

load bwresults_CS1_PM.mat  %load bandwidth data 

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
plotfigs = 0;  % plot if 1, don't plot if 0 
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

for ml = 1:noICs
    
    clearvars -except simlength ml *_ar bw bwm covindex datatrain plotfigs tstep K Fcons hx hy obsfreq *_full LSall CRPSall RMSEall ENSVARall noICs obserrvar ENSMEANERRall  *all_anal
    
    %define the truth here:
    x_obs = x_obs_full(:,:,ml);
    x_true = x_true_full(:,:,ml);
    errx_true = errx_true_full(:,:,ml);

    %--------------------------------------------------------------------------
    %Input parameters:
    n = 1000;  % no. particles or ensemble members 
    x0var = 0; %error variance on x(0)
    novars = K;  %total no. of vars including the obs variable 
    obsind = [3 4 8 9];  %only every second variable observed
    hidind = [1 2 5 6 7];  %hidden states
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
    
    eysample = NaN*ones(novars, n, simlength);  %to store the sampled error distribution
    arvec = NaN*ones(novars, n,simlength);
    arvec(:,:,1) = repmat(errx_true(:,2), 1, n);
    errbnds = [min(datatrain(:,1))-0.05, max(datatrain(:,1))+0.05];
    %--------------------------------------------------------------------------
    %STEP 2. Start the cycling process

    % options=optimset('Display','off');  %this is for fsolve
    % warning off 

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

        fixedres = 1;
        
        for p = 1:novars       
            
            %now sample values:
            clear temperrest 
            
            arvec(p,:,i) = statesA(p,:,i) - tempensstates(p,:,i);
            arvec(p,:,i) = max([arvec(p,:,i); errbnds(1)*ones(1,n)]);
            arvec(p,:,i) = min([arvec(p,:,i); errbnds(2)*ones(1,n)]);
            
            covvaltemp = [statesAnew(p,:,i)', arvec(p,:,i)'];
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
        %STEP 6. Plot results:

        %insert a plotting feature so that we can view the results as we go:

        %first do kernel density estimates for each state:

        %also calculate kernel density estimates so we can plot as we go:
        xpts = repmat([1.5:i-0.5], 2, 1);

        for k = 1:novars
            [ffit(k,:,i), xfit(k,:,i)] = ksdensity(squeeze(ensstates(upstind(k),:,i+1)), 'Weights', wghtall(:,i)); 
       
            [cdff(k,:,i)] = ksdensity(squeeze(ensstates(upstind(k),:,i+1)), xfit(k,:,i), 'function', 'cdf', 'Weights', wghtall(:,i));  %evaluate the cdf also
            normf(k,:,i) = ffit(k,:,i)/max(ffit(k,:,i));

            %going to set values outside the 2.5% and 97.5% range to NaN (ie. to
            %not plot) - the normalising will then only occur for those values that
            %aren't nans
            a = find(or(cdff(k,:,i) < 0.025, cdff(k,:,i) > 0.975));  %ie. outside of this range
            normfplot(k,:,i) = normf(k,:,i);
            if ~isempty(a)
                normfplot(k,a,i) = NaN;
            end       
        end

        if plotfigs == 1

            figure(1)
            clf
            for k = 1:novars
                X = repmat([0.5, xpts(:)', i+0.5], size(xfit,2),1); 
                if i == 1            
                    Y = reshape(repmat(squeeze(xfit(k,:,:))', 2, 1), size(xfit,2), i*2);   
                    C = reshape(repmat(squeeze(normfplot(k,:,:))', 2, 1), size(xfit,2), i*2);
                else
                    Y = reshape(repmat(squeeze(xfit(k,:,:)), 2, 1), size(xfit,2), i*2);   
                    C = reshape(repmat(squeeze(normfplot(k,:,:)), 2, 1), size(xfit,2), i*2);
                end

                subplot(2,5,k)
                h = pcolor(X,Y,C);
                set(h, 'EdgeColor', 'none');
                colorbar        
                hold on 
                plot([1:i], x_true(upstind(k),2:i+1), '.g')
                obsloc = find(obsind == upstind(k));
                if ~isempty(obsloc)
                    hold on
                    plot([1:i], x_obs(obsind(obsloc),2:i+1), '.r')
                end        
                xlabel('time'); ylabel(['x', num2str(k)]);
                ylim auto
            end
        end

        drawnow
    end
   
end


%now calculate the logarithmic score:
zeroval = log(10^-8);  %this is for log score, when the estimated probability is zero.
for k = 1:novars
    [ls(k,:), els(k)] = logscore(squeeze(xfit(k,:,:)), squeeze(ffit(k,:,:)), x_true(k,2:simlength+1), zeroval); %Note, calculating on the truth...
end

save('assimilateresults_CS1_PM.mat')




