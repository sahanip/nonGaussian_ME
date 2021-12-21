%Assimilation CS2, B1

clear all
close all

%--------------------------------------------------------------------------
%first load and extract all the relevant information: 
load evaldata_syn_L96_CS2.mat  %evaluation period data!
tstep = ts;  %can't use h because that's often used as a figure handle.

simlength = 100;  %first 1:simlength part of each IC is run!

x_obs_full = x_obs(:, 1:simlength+1,:);
x_true_full = x_true(:, 1:simlength+1,:);
errx_true_full = errx_true(:, 1:simlength+1,:);

clearvars -except tstep K Fcons hx hy obsfreq *_full noICs obserrvar simlength 

load trainresults_CS2_B1_lambda=2.mat  %B AND PVEC DATA

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
plotfigs = 0;  % plot if 1, don't plot if 0
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX



%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%TUNING PARAMETERS:
%INFLATION:
lambda1 = 1;  %NO INFLATION  %this is equivalent to delta1 
maxinf = 10;  %max val of trace(cov(ensstates)) at which point inflation is turned off - SAME FOR ALL!  
%OPTION3 - GASPARI COHN:
ro_loc = 0; %NO LOCALIZATION  %zero means no localization 
alpha = 0.8;  %tuning parameter for the added model error! 
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

for ml = 1:noICs
    
    clearvars -except simlength ml bmvec driftall Pmmat xt1all lambda1 maxinf ro_loc alpha plotfigs tstep K Fcons hx hy obsfreq *_full LSall CRPSall RMSEall ENSVARall noICs obserrvar ENSMEANERRall *all_anal
    
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
    obserr = obserrvar(obsind,obsind); %obs error statistics (mean and var) known perfectly and assume independence
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

    %LOCALIZATION MATRIX:
    if exist('ro_loc') ~= 0
        if ro_loc == 0
            %no localization 
            CLOC = ones(novars, novars);
        else           
            CLOC = localisation_full_GG(novars,ro_loc);
        end
    elseif isempty(locon)
        CLOC = ones(novars,novars);
    end

    %--------------------------------------------------------------------------
    %STEP 2. Start the cycling process

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

        ensstates(:,:,i+1) = tempensstates(:,:,i+1) - alpha*mvnrnd(repmat(bmvec, 1, n)', repmat(Pmmat, [1, 1, n]))';

        %check physical:
        unph = find(or(ensstates(:,:,i+1) < stbound(:,:,1), ensstates(:,:,i+1) > stbound(:,:,2)));
        while ~isempty(unph)
            enstemp = tempensstates(:,:,i+1) - alpha*mvnrnd(repmat(bmvec, 1, n)', repmat(Pmmat, [1, 1, n]))';
            ensold = ensstates(:,:,i+1);
            ensold(unph) = enstemp(unph);
            ensstates(:,:,i+1) = ensold;
            clear enstemp ensold 
            unph = find(or(ensstates(:,:,i+1) < stbound(:,:,1), ensstates(:,:,i+1) > stbound(:,:,2)));
        end
        %--------------------------------------------------------------------------
        %STEP 4. Undertake DA Update 

        %USE ETKF with inflation and localization:
        Rall(:,:,i) = obserr;
        trcov = sum(diag(cov(ensstates(:,:,i+1)')));
        if trcov > maxinf
            lambdanew = 1;
        else
            lambdanew = lambda1;
        end
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
            %[ffit_anal(k,:,i), xfit_anal(k,:,i)] = ksdensity(squeeze(statesA(upstind(k),:,i+1)), 'Weights', wghtall(:,i+1)); 

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

save(['assimilateresults_CS2_B1_alpha=', num2str(alpha), '.mat'])

