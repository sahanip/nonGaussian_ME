%this creates data for training for case study 1 with increased observation
%error 

clear all
close all

%Lorenz Parameters:
Fcons = 10;
K = 9;  %no. slow vars
eps_g = 1/128;
J = 1/eps_g;  %no. fast vars for each slow variable 
X0tmp = Fcons*rand(K,1); %initial condition (for running on ,not the true IC!)
Y0tmp = rand(J*K,1); %initial condition (for running on ,not the true IC!)
hx = -0.8;
hy = 1; 
ts = 0.0008; 
ts_init = 0.0008;
n=1;
noICs = 30;
trainpoints = 1365;  %total no. of points to be extracted from training period for each separate IC.
evalpoints = 600; %includes the spacing in between successive runs! 
warmuplengthdays = 90;  %in days, then uses 1MTU = 5days to convert. 

obsfreq = 25;  % observe every obsfreq time steps. 
decorrtimestep = 20;  %in terms of no. of obs, no. of obs timesteps needed for autocorr to drop to zero.  Use 20 if obsfreq=25.  This can only be known once a full model run is complete...
obserrmean = zeros(K,1);  %mean zero 
obserrvar = (2*10^-5)*eye(K); %so that model error std is roughly 4 times obs err std %assume independent!
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

%plot errors:
figure
for v = 1:K
    subplot(2,5,v)
    plot(covs_true2(v,1:decorrtimestep:end), errx_true2(v,1:decorrtimestep:end), '.')
    xlabel(['x_', num2str(v), '(t-1)'])
    ylabel(['ex_', num2str(v), '(t)'])
end

%now do transition histogram figures:

figure
counter = 0;
for k = 1:K
    counter = counter+1;
    subplot(2,5,counter)    
    h1 = histogram(errx_true2(k,1:decorrtimestep:end));    
    hold on 
    h2 = histogram(err_obs2(k,1:decorrtimestep:end));     
    title(['x', num2str(k)])
    if k == 1
        legend([h1 h2], 'transition err', 'obs err')
    end
end

counter = 1;
sgtend = NaN*ones(K, size(y_trueall,2), size(y_trueall,3));
for k = 1:K
    sgtend(k,:,:) = (-hx/J)*(sum(y_trueall(counter:counter+J-1,:,:)));
    counter = counter+J;
end

figure(5)
x1_sg = squeeze(x_trueall(1,1:obsfreq*decorrtimestep:end,:));
y1_sg = squeeze(-sgtend(1,1:obsfreq*decorrtimestep:end,:));
plot(x1_sg(:), y1_sg(:), '.')
xlabel('X1'); ylabel('Sub grid tendency')

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%only save selected variables to save space:
clear x_simonestep x_simonetemp x_obsall y1temp z1temp x1temp x_sim h1 h2 errx_trueall x_trueall y_trueall sgtend err_obs2 errx_true2 covs_true2 t_train x1temp  y0_true Y0tmp X0tmp y_true err_obs obstimeind
save('traininputdata_obserr_syn_L96_CS1.mat')

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX








