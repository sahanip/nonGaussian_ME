%script for method B1, case study 1

clear all
close all

%--------------------------------------------------------------------------
%first load and extract all the relevant information:  
load traininputdata_syn_L96_CS1.mat 
tstep = ts;  
x_truefull = x_true;
x_obsfull = x_obs;

clearvars -except Fcons x_truefull x_obsfull J K hx hy c b noICs obserrvar obsfreq tstep

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
plotfigs = 0;  % plot if 1, don't plot if 0
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
%TUNING PARAMETERS:
%INFLATION:
lambda1 = 1.45;   %this is equivalent to delta1 
maxinf = 10;  %max val of trace(cov(ensstates)) at which point inflation is turned off - SAME FOR ALL! 
obserrinflation = 1;   %to inflate particles and accept more 
%LOCALIZATION: GASPARI COHN:
ro_loc = 0;  %zero means no localization 
%$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
    
driftall = []; xt1all = []; driftallv2 = [];

for lm = 1:noICs
    
    lm
    
    clearvars -except counter lm lm1 lm2 rmseanal CLOC ro_locvec ro_loc maxinf driftall driftallv2 xt1all Fcons x_truefull x_obsfull J K hx hy c b noICs lambda1vec lambda1 locon obserrinflation plotfigs lm obserrvar obsfreq tstep funcfilepath

    x_true = x_truefull(:,2:end,lm);
    x_obs = x_obsfull(:,2:end,lm);

    %--------------------------------------------------------------------------
    %Input parameters:
    n = 500;  % no. particles or ensemble members 
    simlength = length(x_obs)-1; %ignore the IC 
    x0var = 0.001; %error variance on x(0)
    novars = K;  %total no. of vars including the obs variable 
    obsind = [3 4 8 9];  %only every second variable observed
    hidind = [1 2 5 6 7];  %hidden states
    H = zeros(length(obsind), novars);  %define observation operator 
    for l = 1:length(obsind)
        H(l, obsind(l)) = 1;
    end
    dtrng = 0.05;
    obserr = obserrvar(obsind,obsind); 
    upstind = [1:K]; %index of internal states to update 
    zeroval = log(10^-8);  %this is for log score, when the estimated probability is zero.
    %--------------------------------------------------------------------------
    %PFMCMC parameters:
    nt = n/2;  %threshold to initiate resampling 
    nx = length(upstind);  %total dimension size 
    hopt = ((4/(nx+2))^(1/(nx+4)))*n^(-1/(nx+4));
    %--------------------------------------------------------------------------
    %Enkf parameters:
    tol = 10^-10;
    %--------------------------------------------------------------------------

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

    %STEP 1. First generate an ensemble of initial conditions, and preallocate
    %necessary matrices:
    x0_true = x_true(:,1);
    x0cov = x0var*eye(novars, novars);
    for k = 1:length(hidind)
        x0cov(hidind(k),hidind(k)) = x0var*50;
    end
    initstates = mvnrnd(repmat(x0_true', n,1), repmat(x0cov, [1, 1, n]))';

    %pre-allocate matrices:
    tempensstates = NaN*ones(novars, n, simlength+1);
    tempensstates(upstind,:,1) = initstates;
    ensstates = tempensstates; statesA = tempensstates; statesAact = tempensstates; statesAnew = statesA;
    Rall = NaN*ones(length(obsind), length(obsind), simlength);

    wghtall = NaN*ones(n, simlength+1);  %to store the weights 
    wghtall(:,1) = 1/n;  %all equally weighted first 

    %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
    %state variable bounds for the PFMCMC - FOLLOWING WILL CHANGE DEPENDING ON THE PROBLEM
    xlb = -25; xub = 25;
    stbound =  xlb*ones(K, n);
    stbound(:,:,2) = xub*ones(K, n);
    stboundeq = NaN*ones(size(stbound));
    stdyncons = [];stdynconseq = [];
    xtbound = stbound(hidind,:,:);  %for the hidden states 
    ybound = stbound(obsind,:,:); %for the observed states
    %XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

    %--------------------------------------------------------------------------
    %STEP 2. Start the cycling process

    % options=optimset('Display','off');  %this is for fsolve
    % warning off 

    for i = 1:simlength 


        cd(funcfilepath)

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
        if isempty(keeppts)
            %unsuccessful run, break
            break
        else

            if length(keeppts) < n
                %redefine initial conditions:          
                statesAnew(:,unph,i) = statesAnew(:, keeppts(ceil(rand([1,length(unph)])*length(keeppts))), i);
                %statesAnew(hidind,:,i) = mvnrnd(statesAnew(hidind,:,i)', repmat(x0cov(hidind,hidind), [1, 1, n]))';    

                %run again:
                sttemp = NaN*ones(K, n, obsfreq+1); sttemp(:,:,1) = statesAnew(:,:,i);  %at the last observed time
                for j = 1:obsfreq
                    sttemp(:,:,j+1) = lorenz96(sttemp(:,:,j),Fcons,K, tstep,n);
                end
                tempensstates(:,:,i+1) = sttemp(:,:,end);
            end
        end
        %--------------------------------------------------------------------------
        %STEP 2. Determine covariates for conditional distribution for each variable:

        ensstates(:,:,i+1) = tempensstates(:,:,i+1);
 
        %--------------------------------------------------------------------------
        %STEP 4. Undertake DA Update 

        %USE ETKF with inflation and localization:
        Rall(:,:,i) = obserr*obserrinflation;
        trcov = sum(diag(cov(ensstates(:,:,i+1)')));
        if trcov > maxinf
            lambdanew = 1;
        else
            lambdanew = lambda1;
        end
        [statesA(:,:,i+1), ~, ~, ensstatesinf(:,:,i+1)] = ETKF_inf_loc(x_obs(obsind,i+1), ensstates(:,:,i+1), Rall(:,:,i), n,H, tol, lambdanew, CLOC);
        wghtall(:,i+1) = 1/n;  %equally weighted for enkf            

        %--------------------------------------------------------------------------
        %STEP 6. Plot results:

        %insert a plotting feature so that we can view the results as we go:

        if plotfigs == 1

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

                subplot(2,4,k)
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

    %save data - only need ensemble mean difference between forecast and
    %analysis:
    driftall = [driftall, squeeze(mean(statesA(:,:,2:end) - ensstates(:,:,2:end),2))];
    %driftallv2 = [driftallv2, squeeze(mean(statesA(:,:,2:end) - ensstatesinf(:,:,2:end),2))];
    xt1all = [xt1all, squeeze(mean(statesA(:,:,1:end-1),2))];    %going to use ensemble mean state also...  
     
end

%now calculate bias and covariance quantities:
bmvec = mean(driftall,2);
Pmmat = cov(driftall');

figure
[x1,x2] = ksdensity(driftall(1,:));  
plot(x2, x1, 'r')
[x5,x6] = ksdensity(driftall(2,:)); 
hold on 
plot(x6, x5, 'g')

figure
plot(xt1all(1,:), driftall(1,:), '.r')%, xt1all(1,:), driftallv2(1,:), '.b')


%save data:
save(['trainresults_CS1_B1_lambda=', num2str(lambda1), '.mat'], 'bmvec', 'Pmmat', 'driftall', 'xt1all')




