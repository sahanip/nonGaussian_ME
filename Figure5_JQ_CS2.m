%this script calcualtes and plots the cost function values (JQ part) for the
%estimated errors of the truth, proposed and 4dvar method using the new Q 

clear all
close all

%CASE STUDY 2
%proposed method:
load alldata_train_CS2_PM.mat

%now change everything so the time indices line up:
errtrue = errx_true;  errestprop = errx_est;  %these are the full ones, not the ones sampled at the decorrelation time. 
xtrue = x_true;
xestprop = NaN*tx_est;  xestprop(:,1:end-1,:) = tx_est(:,2:end,:); 
txestprop = tx_est;
errtruesamp = errx_t_estall;

clearvars -except errtrue *prop basefilepath xtrue windl0 errtruesamp

%method B2
load alldata_train_CS2_B2.mat
errest4dvar = errx_est;
xest4dvar = NaN*tx_est;  xest4dvar(:,1:end-1,:) = tx_est(:,2:end,:); 
txest4dvar = tx_est;
clearvars -except errtrue *prop *4dvar basefilepath obsind K xtrue windl0 errtruesamp Q_est R_true Fcons tstep obsind windl obsfreq simlength x_obs stpt hidind  


Q_true = cov(errtruesamp'); 
setrun = 9;

J4dvar = NaN*ones(1,simlength); Jprop = J4dvar;  Jtrue = Jprop; J4dvar_err = Jprop;  J4dvar_obs = Jprop; Jprop_err  = Jprop;  Jprop_obs = Jprop; Jtrue_err  = Jprop;  Jtrue_obs = Jprop;
RI_true = inv(R_true);
QI_est = inv(Q_true);

for t = stpt:(simlength - windl0 - 30) 
    windl = min(windl0, simlength-t);
    yobswin = x_obs(obsind,t-1:t-1+windl,setrun); 

    xi = xest4dvar(:,t-1,setrun);  
    X = errest4dvar(:,t:t+windl-1,setrun);  %the initial guess/estimate for the errors over the time window 
    [J4dvar(t), ~, Jerrtemp, Jobstemp] = L96testFSOLVEF_window10_2(xi, X, Fcons, K, tstep, yobswin, obsind,windl, obsfreq, RI_true,QI_est); 

    J4dvar_err(t) = sum(Jerrtemp);  J4dvar_obs(t) = sum(Jobstemp);
    clear Jerrtemp Jobstemp 

    xi = xestprop(:,t-1,setrun);
    X = errestprop(:,t:t+windl-1,setrun);  
    [Jprop(t), ~, Jerrtemp, Jobstemp] = L96testFSOLVEF_window10_2(xi, X, Fcons, K, tstep, yobswin, obsind,windl, obsfreq, RI_true,QI_est); 
    Jprop_err(t) = sum(Jerrtemp);  Jprop_obs(t) = sum(Jobstemp);
    clear Jerrtemp Jobstemp 

    xi = xtrue(:,t-1,setrun);
    X = errtrue(:,t:t+windl-1,setrun);
    [Jtrue(t), ~, Jerrtemp, Jobstemp] = L96testFSOLVEF_window10_2(xi, X, Fcons, K, tstep, yobswin, obsind,windl, obsfreq, RI_true,QI_est); 

    Jtrue_err(t) = sum(Jerrtemp);  Jtrue_obs(t) = sum(Jobstemp);
    clear Jerrtemp Jobstemp 
end


figure
h = plot([1:length(Jtrue)], J4dvar_err, 'k', [1:length(Jtrue)], Jprop_err, 'r',[1:length(Jtrue)], Jtrue_err, 'b');
legend('Method B2', 'Proposed','True')
ylabel('{\it{J}}_{\it{Q}}')
xlabel('time step')
xlim([0 500])
set(h, 'LineWidth', 2)
set(gca, 'FontSize', 18,'TickLength',[0.025 0.025], 'LineWidth', 1.8)


