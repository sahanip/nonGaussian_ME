%this script shows the evolution of a solution for the proposed approach.
%intended to serve as a schematic for the proposed approach for case study
%1

clear all
close all

setrun = 7;  
starttime = 2;  

%use the truth as the intiial conditions:
load alldata_train_CS1_PM.mat
initconds = x_obs(:,starttime-1,setrun);
errtruesamp = errx_t_estall;
Q_est = cov(errtruesamp');

obsind = [3 4 8 9];  %only every second variable observed
R_true = obserrvar(obsind,obsind);
RI_true = inv(R_true);
QI_est = inv(Q_est);
clearvars -except setrun starttime initconds basefilepath RI_true QI_est 

%%
%do one iteration of the proposed approach 
load traininputdata_syn_L96_CS1.mat

obsind = [3 4 8 9];  %only every second variable observed
hidind = [1 2 5 6 7];  %hidden states
simlength = size(x_obs,2)-1;%length(t_train);
windl0 = 25;
bwstd =0.5;

stpt = 2;  %starting point for assessment - MUST BE 2, CAN'T CHANGE THIS!
fvalthresh = 0.001;
tstep= ts; clear ts


%--------------------------------------------------------------------------
%Begin cycling:
options=optimset('Display','iter');  %this is for fsolve
%options = optimoptions(@fsolve,'Algorithm','trust-region-dogleg','MaxIterations',1)

warning off 

errx_est = NaN*ones(K, simlength, noICs);
tx_est = errx_est;
fvals_est = NaN*ones(noICs, simlength);

%pre-allocate matrices 
residobs = NaN*ones(length(obsind), simlength); tyobs = residobs; newant = NaN*ones(K, simlength); txsim = newant;  tx = newant; txtrue = newant; txhat = newant; 
errx = zeros*ones(size(newant));  fvals = NaN*ones(1, simlength);
residhid = NaN*ones(length(hidind), simlength); txpsuedoobs = residhid;

newant(:,1) = x_obs(:,1,m);  %going to start with observed initial conditions.
newant(:,starttime-1) = initconds;

m = setrun;
maxiters = 100;

extest = NaN*ones(length(hidind),windl0, maxiters);
extest(:,:,1) = 0;

for t = starttime:starttime
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
    logfunc3 = @(X) L96testFSOLVEF_window2021_withbin(xi, X, Fcons, K, tstep, yobswin, obsind,hidind, windl, obsfreq, bwstd);
    
    [exhidnew, fval] = fsolve(logfunc3, exhid, options); %output gives optimal x(t-2) for yobs   

    [~,~,txtemp, eyall] = L96testFSOLVEF_window2021_withbin(xi, exhidnew, Fcons, K, tstep, yobswin, obsind,hidind, windl, obsfreq, bwstd);  %going to assume true IC known to start off!
       
    %now do the optimization step by step and save results:
    options = optimoptions(@fsolve,'Algorithm','trust-region-dogleg','MaxIterations',1, 'Display', 'off');
    for j = 1:maxiters   %need to choose carefully 
        j
        [extest(:,:,j+1), fvaltest(j+1)] = fsolve(logfunc3, extest(:,:,j), options); %output gives optimal x(t-2) for yobs 
        
        [~,~,~, ~, xbinc(:,j),xbinind(:,j),xbins(j,:)] = L96testFSOLVEF_window2021_withbin(xi, exhidnew, Fcons, K, tstep, yobswin, obsind,hidind, windl, obsfreq, bwstd);
    end
end

%% plotting 

timeplot = [1 2 100];

%in the following, the assumed initial guess on the hidden vars is zero,
%but not necessarily zero for the observed variables.
y1 = errx_true(:,t:t+windl-1);
x1 = x_true(:,t-1:t+windl-2);

ylims = [-0.1 0.1];
xlims = [-6 10];

fs=24

figure
titlelab = {'Initial guess', 'Intermediate', 'Final'};

Posdata = [[0.13 0.229580573951435 0.213405797101449 0.695419426048565]; [0.410797101449275 0.231788079470199 0.213405797101449 0.693211920529801]; [0.691594202898551 0.233995584988962 0.213405797101449 0.691004415011038]];

%tplot = tiledlayout(1,length(timeplot), 'Padding', 'none', 'TileSpacing', 'compact'); 
for l = 1:length(timeplot)
    subplot(1,length(timeplot),l)
   
   % nexttile 
    [~,~,tx1, eyall1,~,~,~,mtest] = L96testFSOLVEF_window2021_withbin(xi, extest(:,:,timeplot(l)), Fcons, K, tstep, yobswin, obsind,hidind, windl, obsfreq, bwstd);
    
   tx2 = tx1(:);
   [~,ind] = sort(tx2);
   condmean = [tx2(ind), mtest(ind)];
    hold on 
    if l == 2
        h1 = plot(x1(:), y1(:), '.r', tx1(:), eyall1(:), '.b', condmean(:,1), condmean(:,2), 'k');
        set(h1(1), 'MarkerSize', 12)
        set(h1(2), 'MarkerSize', 12)
        set(h1(3), 'LineWidth', 2)
    else
        h2 = plot(x1(:), y1(:), '.r', tx1(:), eyall1(:), '.b', condmean(:,1), condmean(:,2), 'k');
        set(h2(1), 'MarkerSize', 12)
        set(h2(2), 'MarkerSize', 12)
        set(h2(3), 'LineWidth', 2)
    end
    ylim(ylims)
    xlim(xlims)
    
    title(titlelab{l})
    if l == 2
        [lh, icons] = legend(h1, {'True', 'Estimates, proposed method', 'Estimated conditional mean'}, 'Location', 'southoutside', 'Orientation', 'horizontal', 'FontSize', 24)
    end
    
    if l ==1
        xlabel('{\it{\bfx}}_{\itj}_-_1[k]')
        ylabel('{\bf{\eta}}_{\itj}[k]')
    end
            
    set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8, 'TitleFontSizeMultiplier', 0.9, 'Position', Posdata(l,:))
    box on 
end

icons = findobj(icons,'Type','line');
% Find lines that use a marker
icons = findobj(icons,'Marker','none','-xor');
% Resize the marker in the legend
set(icons,'MarkerSize',24);

%set(gcf, 'Position',  [100, 100, 800, 200])


annotation(gcf,'textbox',...
    [0.100867724867725 0.844311377245509 0.0214867724867725 0.0718562874251487],...
    'String','a)',...
    'LineStyle','none',...
    'FontSize',fs,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

annotation(gcf,'textbox',...
    [0.383936507936508 0.844311377245507 0.0214867724867725 0.0718562874251487],...
    'String','b)',...
    'LineStyle','none',...
    'FontSize',fs,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');

annotation(gcf,'textbox',...
    [0.667666666666666 0.844311377245505 0.0214867724867726 0.0718562874251487],...
    'String','c)',...
    'LineStyle','none',...
    'FontSize',fs,...
    'FontName','Helvetica Neue',...
    'FitBoxToText','off');




