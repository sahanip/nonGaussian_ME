%execute model error estimation with proposed method for case study 1 and
%increased observation error

clear all
close all

%load training data: 
load traininputdata_obserr_syn_L96_CS1.mat

obsind = [3 4 8 9];  %only every second variable observed
hidind = [1 2 5 6 7];  %hidden states
simlength = size(x_obs,2)-1;%length(t_train);
windl0 = 25;
bwstd =0.5;

stpt = 2;  %starting point for assessment 
fvalthresh = 0.001;
tstep= ts; clear ts

%%
%--------------------------------------------------------------------------
%Begin cycling:
options = optimoptions('fsolve','FunctionTolerance', 10^-3,'Display','off', 'MaxIterations', 200, 'StepTolerance', 10^-3)
warning off 

errx_est = NaN*ones(K, simlength, noICs);
tx_est = errx_est;
fvals_est = NaN*ones(noICs, simlength);

for m = 1:noICs
    
    %pre-allocate matrices 
    residobs = NaN*ones(length(obsind), simlength); tyobs = residobs; newant = NaN*ones(K, simlength); txsim = newant;  tx = newant; txtrue = newant; txhat = newant; 
    errx = zeros*ones(size(newant));  fvals = NaN*ones(1, simlength);
    residhid = NaN*ones(length(hidind), simlength); txpsuedoobs = residhid;

    newant(:,1) = x_obs(:,1,m);  %going to start with observed initial conditions.

    for t = stpt:(simlength - windl0 - 1)
        [m,t]
        clear Sfinal Xhat_t2 Xstar_t2 Xhat_t1 Xstar_t2vec fvalvec Xhat_t2vec ind 

        %--------------------------------------------------------------------------
        %STEP1. Define inputs for optimising initial condition:

        %define inputs:
        windl = min(windl0, simlength-t);
        yobswin = x_obs(obsind,t-1:t-1+windl,m);  %observations over time window 
        xi = newant(:,t-1); %best estimate of initial condition 
        exhid = errx(hidind,t:t+windl-1);  %initial guess of errors for hidden states
        %--------------------------------------------------------------------------
        %STEP2. Estimate the true initial condition for time t:

        %Specify the function to minimise:
        logfunc3 = @(X) L96testFSOLVEF_window2021(xi, X, Fcons, K, tstep, yobswin, obsind,hidind, windl, obsfreq, bwstd);  %going to assume true IC known to start off!
       
        [exhidnew, fval] = fsolve(logfunc3, exhid, options); %output gives optimal x(t-2) for yobs   

        %--------------------------------------------------------------------------
        %STEP 3: Calculate x(t) using optimised new errors:
        
        xsmtemp = NaN*ones(K, obsfreq+1); xsmtemp(:,1) = xi;
        for v = 1:obsfreq
            xsmtemp(:,v+1) =  lorenz96(xsmtemp(:,v),Fcons,K, tstep,1); 
        end
        Xsim_t = xsmtemp(:,end);
        exobsnew = x_obs(obsind,t,m) - Xsim_t(obsind);
        Xhat_t = NaN*ones(K,1);
        Xhat_t(obsind) = x_obs(obsind,t,m);
        Xhat_t(hidind) = Xsim_t(hidind) + exhidnew(:,1);

        %--------------------------------------------------------------------------    
        %STEP 4: Calculate Error in all variables:

        errx(obsind,t) = exobsnew;  
        %input error for hidden states at future times, even though this will be re-written (still used as starting point):
        errx(hidind,t:t+windl-1) = exhidnew;   %estimated transition error 

        %Output the corresponding states for covariate purposes:
        tx(:,t) = xi;  %conditioning covariate for that residual. (x(t-1))
        txtrue(:,t) = x_true(:,t-1,m);  %true conditioning covaraite for that residual.        
        txsim(:,t) = Xsim_t; % one step simulation using best guess for initial condition 
        txhat(:,t) = Xhat_t;  %best estimate at time t (observed variables should equal obs)
        tyobs(:,t) = x_obs(obsind,t,m);   %
        fvals(1,t) = fval;

        %insert the new antecedent states for use in next time step:
        newant(:,t) = Xhat_t; 
        
        %diagnostic purposes:
        [errx(:,t), errx_true(:,t,m)]

    end
    
    %now save data that's required:
    errx_est(:,:,m) = errx;
    tx_est(:,:,m) = tx;
    fvals_est(m,:) = fvals;
    
end

sampselec = [stpt:decorrtimestep:(simlength - windl0 - 1)];  %first 74 values, and sample at decorrelation time step.

%--------------------------------------------------------------------------
%Now collate data for examination:
errx_estall = NaN*ones(K, (length(sampselec))*size(errx_est,3)); tx_estall = errx_estall;
errx_t_estall = errx_estall; tx_t_estall = errx_estall;  obserr_t_estall = errx_estall;  

for v = 1:K
    x1tmp = squeeze(errx_est(v,sampselec,:)); 
    x2tmp = squeeze(tx_est(v,sampselec,:));
    errx_estall(v,:) = x1tmp(:);
    tx_estall(v,:) = x2tmp(:);
    
    x1_t_tmp = squeeze(errx_true(v,sampselec,:));
    x2_t_tmp = squeeze(covs_true(v,sampselec,:));
    errx_t_estall(v,:) = x1_t_tmp(:);
    tx_t_estall(v,:) = x2_t_tmp(:);
    
    x3_t_tmp = squeeze(x_obs(v,sampselec,:) - x_true(v,sampselec,:));
    obserr_t_estall(v,:) = x3_t_tmp(:);
end

%Create data for CKDE:
varsel = 9;

alldata = [errx_estall(varsel,:)', tx_estall(varsel,:)'];
alldata_true = [errx_t_estall(varsel,:)', tx_t_estall(varsel,:)'];
   
%NOTE - NEED TO REMOVE NAN, HENCE END-1
if ~isempty(find(isnan(alldata)))
    disp('NANS')
end

save('alldata_train_obserr_CS1_PM.mat')
save('trainresults_obserr_CS1_PM_forR.mat', 'alldata')
save('trainresults_obserr_CS1_PM_forR_true.mat', 'alldata_true')

%save figures:
figure(1)
saveas(gcf, 'Est_true_errx_xt1.fig')
figure(3)
saveas(gcf, 'Diff_est_xt1.fig')




