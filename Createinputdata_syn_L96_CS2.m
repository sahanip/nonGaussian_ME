%this script constructs a synthetic error model and set of synthetic data
%on which to do the model error estimation, including truth, synthetic obs
%and training data for CS2 

clear all
close all

%Lorenz Parameters:
Fcons = 14;
K = 9;  %no. slow vars
eps_g = 0.7;
J = 20;  %no. fast vars for each slow variable 
hx = -2;
hy = 1; 

X0tmp = Fcons*rand(K,1); %initial condition (for running on ,not the true IC!)
Y0tmp = rand(J*K,1); %initial condition (for running on ,not the true IC!)
ts = 0.0008; 
ts_init = 0.0008;
n=1;
noICs = 30;
trainpoints = 561;  %total no. of points to be extracted from training period for each separate IC.
evalpoints = 600; 
warmuplengthdays = 90;  %in days, then uses 1MTU = 5days to convert. 

obsfreq = 50;  % observe every obsfreq time steps. 
decorrtimestep = 8;  %in terms of no. of obs, no. of obs timesteps needed for autocorr to drop to zero.  This can only be known once a full model run is complete...
obserrmean = zeros(K,1);  %mean zero 
obserrvar = 0.0000001*eye(K); %0.00000001*eye(K);  %assume independent!
stbound = [-15*ones(K,1), 15*ones(K,1)];  %going to force the states to be bounded so as to prevent instability of dynamical system.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create data for training period:

warmuplength = (90/5)*(1/ts_init);

%STEP 1: Single model run to extract ICs:
runtime = (obsfreq*trainpoints);
Xtemp = X0tmp;  Ytemp = Y0tmp;

%first run on to attractor:
for t = 1:warmuplength
    %warm up period:    
    [Xtemp, Ytemp] = lorenz96_eps_fin(Xtemp,Ytemp, Fcons,K,J, ts_init, eps_g, hx, hy);
end

%now run and save different ICs:
x0_true = NaN*ones(K, noICs); y0_true = NaN*ones(J*K, noICs); 
x0_true(:,1) = Xtemp;  y0_true(:,1) = Ytemp;

for k = 1:noICs
    k
    for t = 1:runtime        
        [Xtemp, Ytemp] = lorenz96_eps_fin(Xtemp,Ytemp, Fcons,K,J, ts_init, eps_g, hx, hy);
    end
    if k == noICs
        x0_true_eval = Xtemp; y0_true_eval = Ytemp;
    else
        x0_true(:,k+1) = Xtemp;  y0_true(:,k+1) = Ytemp;
    end    
end

clear Xtemp Ytemp 

%STEP 2: Generate true and model runs for each IC:
%Initialise:
t_train = 0:ts:(obsfreq*trainpoints)*ts;
x_sim = NaN*ones(K, length(t_train), noICs); x_sim(:,1,:) = x0_true;  x_obsall = x_sim; x_simonestep = x_sim; errx_trueall = x_sim; 
x_trueall = NaN*ones(K, length(t_train), noICs); x_trueall(:,1,:) = x0_true;
y_trueall = NaN*ones(J*K, length(t_train), noICs); y_trueall(:,1,:) = y0_true;

for m = 1:noICs
%    counter = 1;
    m
    for i = 1:length(t_train)-1
        %i        
        %assumed model
        x_sim(:,i+1,m) = lorenz96(x_sim(:,i,m),Fcons,K, ts,n);

        %one step version (i.e. initial condition correct, transition incorrect):
        x_simonestep(:,i+1,m) = lorenz96(x_trueall(:,i,m),Fcons,K, ts,n);
        
        [x_trueall(:,i+1,m), y_trueall(:,i+1,m)] = lorenz96_eps_fin(x_trueall(:,i,m),y_trueall(:,i,m), Fcons,K,J, ts, eps_g, hx, hy);

        errx_trueall(:,i+1,m) = x_trueall(:,i+1,m) - x_simonestep(:,i+1,m);

        %now create synthetic obs:
        x_obsall(:,i+1,m) = x_trueall(:, i+1,m) + mvnrnd(obserrmean', obserrvar)';    

    end
end


%now determine the transition errors that would be used in the training
%period:
obstimeind = [1:obsfreq:length(t_train)];
x_obs = x_obsall(:,obstimeind,:);
err_obs = x_trueall(:, obstimeind,:) - x_obs;

%now generate true errors at observation time step:
errx_true = NaN*ones(K, length(obstimeind), noICs); covs_true = errx_true; x_true = errx_true; y_true = NaN*ones(K*J, length(obstimeind), noICs);
for m = 1:noICs
    for i = 2:length(obstimeind)
        x_simonetemp = NaN*ones(K, obsfreq+1); x_simonetemp(:,1) = x_trueall(:, obstimeind(i-1),m);
        
        for v = 1:obsfreq            
            x_simonetemp(:,v+1) = lorenz96(x_simonetemp(:,v),Fcons,K, ts,n);
        end
        
        errx_true(:,i,m) = x_trueall(:,obstimeind(i),m) - x_simonetemp(:,end);
        covs_true(:,i,m) = x_trueall(:,obstimeind(i-1),m);
        x_true(:,i,m) = x_trueall(:, obstimeind(i), m);
        y_true(:,i,m) = y_trueall(:, obstimeind(i), m);
    end
end


%now reassign into single vector:
errx_true2 = NaN*ones(K, (size(errx_true,2)-1)*size(errx_true,3)); covs_true2 = errx_true2; err_obs2 = errx_true2;
for v = 1:K
    clear x1temp y1temp
    x1temp  = squeeze(errx_true(v,2:end,:));
    y1temp  = squeeze(covs_true(v,2:end,:));
    z1temp  = squeeze(err_obs(v,2:end,:));
    errx_true2(v,:) = x1temp(:);
    covs_true2(v,:) = y1temp(:);
    err_obs2(v,:) = z1temp(:);
end

%to analyse sub-grid tendencies 
counter = 1;
sgtend = NaN*ones(K, size(y_trueall,2), size(y_trueall,3));
for k = 1:K
    sgtend(k,:,:) = (-hx/J)*(sum(y_trueall(counter:counter+J-1,:,:)));
    counter = counter+J;
end
x1_sg = squeeze(x_trueall(1,1:obsfreq*decorrtimestep:end,:));
y1_sg = squeeze(-sgtend(1,1:obsfreq*decorrtimestep:end,:));


%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%only save selected variables to save space:
clear x_simonestep x_simonetemp x_obsall y1temp z1temp x1temp x_sim h1 h2 errx_trueall x_trueall y_trueall sgtend err_obs2 errx_true2 covs_true2 t_train x1temp  y0_true Y0tmp X0tmp y_true err_obs obstimeind
save('traininputdata_syn_L96_CS2.mat')
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


%% 

close all
clearvars -except warmuplength decorrtimestep Fcons J K hx hy eps_g c b d n ts testno obsfreq obserrmean obserrvar stbound  cbmult ts_init noICs ICspacing evalpoints  expname 

%initial conditions:
X0tmp = Fcons*rand(K,1); %initial condition (for running on ,not the true IC!)
Y0tmp = rand(J*K,1); %initial condition (for running on ,not the true IC!)

%get onto attractor:
Xtemp = X0tmp;  Ytemp = Y0tmp;
for t = 1:warmuplength
    %warm up period:    
    [Xtemp, Ytemp] = lorenz96_eps_fin(Xtemp,Ytemp, Fcons,K,J, ts_init, eps_g, hx, hy);
end

x0_true = Xtemp;
y0_true = Ytemp;
clear Xtemp Ytemp 


%STEP 2: Generate true and model runs for each IC:
%Initialise:
t_eval = 0:ts:(obsfreq*evalpoints)*ts;
x_sim = NaN*ones(K, length(t_eval), noICs); x_sim(:,1,1) = x0_true;  x_obsall = x_sim; x_simonestep = x_sim; errx_trueall = x_sim; 
x_trueall = NaN*ones(K, length(t_eval), noICs); x_trueall(:,1,1) = x0_true;
y_trueall = NaN*ones(J*K, length(t_eval), noICs); y_trueall(:,1,1) = y0_true;

for m = 1:noICs
%    counter = 1;
    m
    for i = 1:length(t_eval)-1
        %i
        if and(m > 1, i == 1)
            i1 = length(t_eval);   
            m1 = m-1;
        else
            i1 = i;
            m1 = m;
        end
            
        %assumed model
        x_sim(:,i+1,m) = lorenz96(x_sim(:,i1,m1),Fcons,K, ts,n);

        %one step version (i.e. initial condition correct, transition incorrect):
        x_simonestep(:,i+1,m) = lorenz96(x_trueall(:,i1,m1),Fcons,K, ts,n);
       
        [x_trueall(:,i+1,m), y_trueall(:,i+1,m)] = lorenz96_eps_fin(x_trueall(:,i1,m1),y_trueall(:,i1,m1), Fcons,K,J, ts, eps_g, hx, hy);

        errx_trueall(:,i+1,m) = x_trueall(:,i+1,m) - x_simonestep(:,i+1,m);

        %now create synthetic obs:
        x_obsall(:,i+1,m) = x_trueall(:, i+1,m) + mvnrnd(obserrmean', obserrvar)';    

    end
end

%now determine the transition errors that would be used in the training
%period:
obstimeind = [1:obsfreq:length(t_eval)];
%x_obs = x_obsall(:,obstimeind,:);

%now generate true errors at observation time step:
errx_true = NaN*ones(K, length(obstimeind)-1, noICs); covs_true = errx_true; x_true = errx_true; y_true = NaN*ones(K*J, length(obstimeind)-1, noICs); x_obs = errx_true;
for m = 1:noICs
    for i = 2:length(obstimeind)
        x_simonetemp = NaN*ones(K, obsfreq+1); x_simonetemp(:,1) = x_trueall(:, obstimeind(i-1),m);
        
        for v = 1:obsfreq            
            x_simonetemp(:,v+1) = lorenz96(x_simonetemp(:,v),Fcons,K, ts,n);
        end
        
        errx_true(:,i-1,m) = x_trueall(:,obstimeind(i),m) - x_simonetemp(:,end);
        covs_true(:,i-1,m) = x_trueall(:,obstimeind(i-1),m);
        x_true(:,i-1,m) = x_trueall(:, obstimeind(i), m);
        y_true(:,i-1,m) = y_trueall(:, obstimeind(i), m);
        x_obs(:,i-1,m) = x_obsall(:, obstimeind(i), m);
    end
end

err_obs = x_true- x_obs;

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
clear x_trueall y_trueall errx_trueall x_obsall t_eval x_sim x_simonestep Y0tmp y_true 
save('evaldata_syn_L96_CS2.mat')
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX













