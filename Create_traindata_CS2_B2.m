%train for B2 CS2  

clear all
close all

%load the Q estimate which will be used in the cost function:
load trainresults_CS2_B1_lambda=2.mat
Q_est = Pmmat;
clearvars -except Q_est

load traininputdata_syn_L96_CS2.mat

testno = 'Test19'

obsind = [1 2 5 6];  %only every second variable observed
hidind = [3 4 7 8 9];  %hidden states
simlength = size(x_obs,2)-1;%length(t_train);
windl0 = 25;

R_true = obserrvar(obsind,obsind);

RI_true = inv(R_true);
QI_est = inv(Q_est);

stpt = 2;  %starting point for assessment 
fvalthresh = 0.001;
tstep= ts; clear ts

%%
%--------------------------------------------------------------------------
%Begin cycling:
options=optimset('Display','off');  %this is for fsolve
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
        exguess = errx(:,t:t+windl-1);  %initial guess of errors for hidden states
        %--------------------------------------------------------------------------
        %STEP2. Estimate the true initial condition for time t:

        %Specify the function to minimise:
        logfunc3 = @(X) L96testFSOLVEF_window10(xi, X, Fcons, K, tstep, yobswin, obsind,windl, obsfreq, RI_true,QI_est);  %going to assume true IC known to start off!

        [exguessnew, fval] = fsolve(logfunc3, exguess, options); %output gives optimal x(t-2) for yobs   

        %--------------------------------------------------------------------------
        %STEP 3: Calculate x(t) using optimised new errors:
        
        xsmtemp = NaN*ones(K, obsfreq+1); xsmtemp(:,1) = xi;
        for v = 1:obsfreq
            xsmtemp(:,v+1) =  lorenz96(xsmtemp(:,v),Fcons,K, tstep,1); 
        end
        Xsim_t = xsmtemp(:,end);
        Xhat_t = Xsim_t + exguessnew(:,1);

        %--------------------------------------------------------------------------    
        %STEP 4: Calculate Error in all variables:

        errx(:,t:t+windl-1) = exguessnew;

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

norows = 2 
nocols = 5

sampselec = [stpt+1:(simlength - windl0 - 1)]; 
%specify no. points to select:
nopts = 6000;

%--------------------------------------------------------------------------
%Now collate data for examination:
errx_estall = NaN*ones(K, (length(sampselec))*size(errx_est,3)); tx_estall = errx_estall; errxt1_estall = errx_estall;
errx_t_estall = errx_estall; tx_t_estall = errx_estall;  errxt1_t_estall = errx_estall;

for v = 1:K
    x1tmp = squeeze(errx_est(v,sampselec,:));  %TO AVOID NANS FROM LAST STEP!
    x2tmp = squeeze(tx_est(v,sampselec,:));
    x3tmp = squeeze(errx_est(v,sampselec-1,:)); 
    errx_estall(v,:) = x1tmp(:);
    tx_estall(v,:) = x2tmp(:);
    errxt1_estall(v,:) = x3tmp(:);
    
    x1_t_tmp = squeeze(errx_true(v,sampselec,:));
    x2_t_tmp = squeeze(covs_true(v,sampselec,:));
    x3_t_tmp = squeeze(errx_true(v,sampselec-1,:));
    errx_t_estall(v,:) = x1_t_tmp(:);
    tx_t_estall(v,:) = x2_t_tmp(:);
    errxt1_t_estall(v,:) = x3_t_tmp(:);
    
end

%Create data for CKDE:
varsel = 2;

%now first calclate approx. density of covariates:
ballsize = [1,1,0.01]; %for xk, xk-1,exk directions.
fullcovdata = [tx_estall(varsel,:)', tx_estall(varsel-1,:)', errxt1_estall(varsel,:)'];
for k = 1:size(fullcovdata,1)
    disttemp = abs(fullcovdata  - repmat(fullcovdata(k,:), size(fullcovdata,1),1));
    ptdens(k) = length(find(((disttemp(:,1) < ballsize(1)) + (disttemp(:,2) < ballsize(2)) + (disttemp(:,3) < ballsize(3))) == 3));
end


%first select points with the desired density of points:
densthresh = 20; %no points within ball for it to be considered acceptable:
min_e = -0.38; max_e = 0.15;
pts1 = find(ptdens > densthresh);
keeppts = pts1(randsample([1:length(pts1)], nopts));

%3covs:
alldata = [errx_estall(varsel,keeppts)', tx_estall(varsel,keeppts)', tx_estall(varsel-1,keeppts)', errxt1_estall(varsel,keeppts)'];


%NOTE - NEED TO REMOVE NAN, HENCE END-1
if ~isempty(find(isnan(alldata)))
    disp('NANS')
end

clear datanew ebins eptsreq ecount ebinnew emax 

save('alldata_train_CS2_B2.mat')
save('trainresults_CS2_B2_forR.mat', 'alldata')

