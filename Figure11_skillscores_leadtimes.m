%Figure: skill scores for different lead times for all 3 methods
%method for Case study 1.  Plots CRPS, LS, RMSE.

%separates into obs and hidden, and also calculates variability by
%considering the mean for each simulation run.

clear all
close all

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%USER INPUTS:

%CASE STUDY 1:
%following is to load obsind and hidind:
load assimilateresults_CS1_PM.mat
x_true_cs1 = x_true;
clearvars -except *_cs1 

%load proposed 
load forecastresults_CS1_PM.mat
CRPSall1_cs1 = CRPSall; ENSVARall1_cs1 = ENSVARall; LSall1_cs1 = LSall; ENSMEANERRall1_cs1 = ENSMEANERRall;
clearvars -except *_cs1

%now method b1 
load forecastresults_CS1_B1_alpha=0.8.mat
CRPSall2_cs1 = CRPSall; ENSVARall2_cs1 = ENSVARall; LSall2_cs1 = LSall; ENSMEANERRall2_cs1 = ENSMEANERRall;
clearvars -except *_cs1

%Now themehtod b2 
load forecastresults_CS1_B2.mat
CRPSall3_cs1 = CRPSall; ENSVARall3_cs1 = ENSVARall; LSall3_cs1 = LSall; ENSMEANERRall3_cs1 = ENSMEANERRall;
maxleadtime_cs1 = maxleadtime;
clearvars -except *_cs1  

%-------------------
%CASE STUDY 2:
%following is to load obsind and hidind:
load assimilateresults_CS2_PM.mat
x_true_cs2 = x_true;
clearvars -except *_cs1 *_cs2

%load proposed approach 
load forecastresults_CS2_PM.mat
CRPSall1_cs2 = CRPSall; ENSVARall1_cs2 = ENSVARall; LSall1_cs2 = LSall; ENSMEANERRall1_cs2 = ENSMEANERRall;
clearvars -except *_cs1 *_cs2

%now method b1 
load forecastresults_CS2_B1_alpha=0.8.mat
CRPSall2_cs2 = CRPSall; ENSVARall2_cs2 = ENSVARall; LSall2_cs2 = LSall; ENSMEANERRall2_cs2 = ENSMEANERRall;
clearvars -except *_cs1 *_cs2

%Now method b2 
%load my results scores:
load forecastresults_CS2_B2.mat
CRPSall3_cs2 = CRPSall; ENSVARall3_cs2 = ENSVARall; LSall3_cs2 = LSall; ENSMEANERRall3_cs2 = ENSMEANERRall;
maxleadtime_cs2 = maxleadtime;
clearvars -except *_cs1 *_cs2 

legnames = {'B1 Method', 'B2 Method'}; %order: 2,3
mtu_cs1 = 0.02;  %mtu of spacing of observations 
mtu_cs2 = 0.04;  %mtu of spacing of observations 

perfect_RMSE = 0;
perfect_CRPS = 0;
perfect_lsbuff = 0.05;  %amount less than the minimal value observed.

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

%% CASE STUDY 1:
noICs = length(CRPSall1_cs1);

%first aggregate the data and calculate skill scores:

CRPS1_cs1 = cell(1,maxleadtime_cs1); CRPS2_cs1 = CRPS1_cs1; CRPS3_cs1 = CRPS1_cs1; CRPSS2_cs1 = CRPS1_cs1;  CRPSS3_cs1 = CRPS1_cs1;
LS1_cs1 = cell(1,maxleadtime_cs1); LS2_cs1 = LS1_cs1; LS3_cs1 = LS1_cs1; LSS2_cs1 = LS1_cs1;  LSS3_cs1 = LS1_cs1;
ENSMEANERR1_cs1 = cell(1,maxleadtime_cs1); ENSMEANERR2_cs1 = ENSMEANERR1_cs1; ENSMEANERR3_cs1 = ENSMEANERR1_cs1; 
RMSE1_cs1 = cell(1,maxleadtime_cs1); RMSE2_cs1 = RMSE1_cs1; RMSE3_cs1 = RMSE1_cs1; RMSES2_cs1 = RMSE1_cs1;  RMSES3_cs1 = RMSE1_cs1;

CRPSS2mean_cs1 = NaN*ones(noICs, maxleadtime_cs1); CRPSS3mean_cs1 = CRPSS2mean_cs1;
LSS2mean_cs1 = NaN*ones(noICs, maxleadtime_cs1); LSS3mean_cs1 = LSS2mean_cs1;
RMSES2mean_cs1 = NaN*ones(noICs, maxleadtime_cs1); RMSES3mean_cs1 = RMSES2mean_cs1;

%first calculate the min LS:
minls1 = NaN*ones(1, noICs); minls2 = NaN*ones(1, noICs); minls3 = NaN*ones(1, noICs);
for m = 1:noICs
    minls1(m) = min(min(min(LSall1_cs1{m})));
    minls2(m) = min(min(min(LSall2_cs1{m})));
    minls3(m) = min(min(min(LSall3_cs1{m})));
end

perfect_LS_cs1 = min([minls1, minls2, minls3]) - perfect_lsbuff;
    

for k = 1:maxleadtime_cs1
    clear crps*temp ls*temp ens*temp rmse*tmp
    for m = 1:noICs
        crps1temp(:,:,m) = squeeze(CRPSall1_cs1{m}(:,:,k)); crps2temp(:,:,m) = squeeze(CRPSall2_cs1{m}(:,:,k)); crps3temp(:,:,m) = squeeze(CRPSall3_cs1{m}(:,:,k));
        ls1temp(:,:,m) = squeeze(LSall1_cs1{m}(:,:,k)); ls2temp(:,:,m) = squeeze(LSall2_cs1{m}(:,:,k)); ls3temp(:,:,m) = squeeze(LSall3_cs1{m}(:,:,k));
        ensmeanerr1temp(:,:,m) = squeeze(ENSMEANERRall1_cs1{m}(:,:,k)); ensmeanerr2temp(:,:,m) = squeeze(ENSMEANERRall2_cs1{m}(:,:,k)); ensmeanerr3temp(:,:,m) = squeeze(ENSMEANERRall3_cs1{m}(:,:,k));
        
        %now calculate rmse:
        rmse1temp(m) = sqrt(nanmean(nanmean((ensmeanerr1temp(:,:,m)).^2)));
        rmse2temp(m) = sqrt(nanmean(nanmean((ensmeanerr2temp(:,:,m)).^2)));
        rmse3temp(m) = sqrt(nanmean(nanmean((ensmeanerr3temp(:,:,m)).^2)));
    end
    CRPS1_cs1{k} = crps1temp; CRPS2_cs1{k} = crps2temp; CRPS3_cs1{k} = crps3temp; 
    LS1_cs1{k} = ls1temp; LS2_cs1{k} = ls2temp; LS3_cs1{k} = ls3temp;     
    ENSMEANERR1_cs1{k} = ensmeanerr1temp; ENSMEANERR2_cs1{k} = ensmeanerr2temp; ENSMEANERR3_cs1{k} = ensmeanerr3temp; 
    RMSE1_cs1{k} = rmse1temp; RMSE2_cs1{k} = rmse2temp; RMSE3_cs1{k} = rmse3temp;
    
    %calculate skill score:
    CRPSS2_cs1{k} = (crps1temp - crps2temp)./(perfect_CRPS - crps2temp); CRPSS3_cs1{k} = (crps1temp - crps3temp)./(perfect_CRPS - crps3temp);    
    LSS2_cs1{k} = (ls1temp - ls2temp)./(perfect_LS_cs1 - ls2temp); LSS3_cs1{k} = (ls1temp - ls3temp)./(perfect_LS_cs1 - ls3temp);
    RMSES2_cs1{k} = (rmse1temp - rmse2temp)./(perfect_RMSE - rmse2temp); RMSES3_cs1{k} = (rmse1temp - rmse3temp)./(perfect_RMSE - rmse3temp);
        
    %now calculate stats - mean over all variables and all time:
    for m = 1:noICs
        x1 = CRPSS2_cs1{k}(:,:,m);
        CRPSS2mean_cs1(m,k) = nanmean(x1(:));      
        x1 = CRPSS3_cs1{k}(:,:,m);
        CRPSS3mean_cs1(m,k) = nanmean(x1(:));   
        
        x1 = LSS2_cs1{k}(:,:,m);
        LSS2mean_cs1(m,k) = nanmean(x1(:));      
        x1 = LSS3_cs1{k}(:,:,m);
        LSS3mean_cs1(m,k) = nanmean(x1(:)); 
        
        RMSES2mean_cs1(m,k) = (RMSES2_cs1{k}(m));      
        RMSES3mean_cs1(m,k) = (RMSES3_cs1{k}(m)); 
    end        
    
end

%% CASE STUDY 2:

noICs = length(CRPSall1_cs2);

%first aggregate the data and calculate skill scores:

CRPS1_cs2 = cell(1,maxleadtime_cs2); CRPS2_cs2 = CRPS1_cs2; CRPS3_cs2 = CRPS1_cs2; CRPSS2_cs2 = CRPS1_cs2;  CRPSS3_cs2 = CRPS1_cs2;
LS1_cs2 = cell(1,maxleadtime_cs2); LS2_cs2 = LS1_cs2; LS3_cs2 = LS1_cs2; LSS2_cs2 = LS1_cs2;  LSS3_cs2 = LS1_cs2;
ENSMEANERR1_cs2 = cell(1,maxleadtime_cs2); ENSMEANERR2_cs2 = ENSMEANERR1_cs2; ENSMEANERR3_cs2 = ENSMEANERR1_cs2; 
RMSE1_cs2 = cell(1,maxleadtime_cs2); RMSE2_cs2 = RMSE1_cs2; RMSE3_cs2 = RMSE1_cs2; RMSES2_cs2 = RMSE1_cs2;  RMSES3_cs2 = RMSE1_cs2;

CRPSS2mean_cs2 = NaN*ones(noICs, maxleadtime_cs2); CRPSS3mean_cs2 = CRPSS2mean_cs2;
LSS2mean_cs2 = NaN*ones(noICs, maxleadtime_cs2); LSS3mean_cs2 = LSS2mean_cs2;
RMSES2mean_cs2 = NaN*ones(noICs, maxleadtime_cs2); RMSES3mean_cs2 = RMSES2mean_cs2;

%first calculate the min LS:
minls1 = NaN*ones(1, noICs); minls2 = NaN*ones(1, noICs); minls3 = NaN*ones(1, noICs);
for m = 1:noICs
    minls1(m) = min(min(min(LSall1_cs2{m})));
    minls2(m) = min(min(min(LSall2_cs2{m})));
    minls3(m) = min(min(min(LSall3_cs2{m})));
end

perfect_LS_cs2 = min([minls1, minls2, minls3]) - perfect_lsbuff;

for k = 1:maxleadtime_cs2
    clear crps*temp ls*temp ens*temp rmse*tmp
    for m = 1:noICs
        crps1temp(:,:,m) = squeeze(CRPSall1_cs2{m}(:,:,k)); crps2temp(:,:,m) = squeeze(CRPSall2_cs2{m}(:,:,k)); crps3temp(:,:,m) = squeeze(CRPSall3_cs2{m}(:,:,k));
        ls1temp(:,:,m) = squeeze(LSall1_cs2{m}(:,:,k)); ls2temp(:,:,m) = squeeze(LSall2_cs2{m}(:,:,k)); ls3temp(:,:,m) = squeeze(LSall3_cs2{m}(:,:,k));
        ensmeanerr1temp(:,:,m) = squeeze(ENSMEANERRall1_cs2{m}(:,:,k)); ensmeanerr2temp(:,:,m) = squeeze(ENSMEANERRall2_cs2{m}(:,:,k)); ensmeanerr3temp(:,:,m) = squeeze(ENSMEANERRall3_cs2{m}(:,:,k));
        
        %now calculate rmse:
        rmse1temp(m) = sqrt(nanmean(nanmean((ensmeanerr1temp(:,:,m)).^2)));
        rmse2temp(m) = sqrt(nanmean(nanmean((ensmeanerr2temp(:,:,m)).^2)));
        rmse3temp(m) = sqrt(nanmean(nanmean((ensmeanerr3temp(:,:,m)).^2)));
    end
    CRPS1_cs2{k} = crps1temp; CRPS2_cs2{k} = crps2temp; CRPS3_cs2{k} = crps3temp; 
    LS1_cs2{k} = ls1temp; LS2_cs2{k} = ls2temp; LS3_cs2{k} = ls3temp; 
    ENSMEANERR1_cs2{k} = ensmeanerr1temp; ENSMEANERR2_cs2{k} = ensmeanerr2temp; ENSMEANERR3_cs2{k} = ensmeanerr3temp; 
    RMSE1_cs2{k} = rmse1temp; RMSE2_cs2{k} = rmse2temp; RMSE3_cs2{k} = rmse3temp;
    
    %calculate skill score:
    CRPSS2_cs2{k} = (crps1temp - crps2temp)./(perfect_CRPS - crps2temp); CRPSS3_cs2{k} = (crps1temp - crps3temp)./(perfect_CRPS - crps3temp);    
    LSS2_cs2{k} = (ls1temp - ls2temp)./(perfect_LS_cs2 - ls2temp); LSS3_cs2{k} = (ls1temp - ls3temp)./(perfect_LS_cs2 - ls3temp);
    RMSES2_cs2{k} = (rmse1temp - rmse2temp)./(perfect_RMSE - rmse2temp); RMSES3_cs2{k} = (rmse1temp - rmse3temp)./(perfect_RMSE - rmse3temp);
        
    %now calculate stats - mean over all variables and all time:
    for m = 1:noICs
        x1 = CRPSS2_cs2{k}(:,:,m);
        CRPSS2mean_cs2(m,k) = nanmean(x1(:));      
        x1 = CRPSS3_cs2{k}(:,:,m);
        CRPSS3mean_cs2(m,k) = nanmean(x1(:));   
        
        x1 = LSS2_cs2{k}(:,:,m);
        LSS2mean_cs2(m,k) = nanmean(x1(:));      
        x1 = LSS3_cs2{k}(:,:,m);
        LSS3mean_cs2(m,k) = nanmean(x1(:)); 
        
        RMSES2mean_cs2(m,k) = (RMSES2_cs2{k}(m));      
        RMSES3mean_cs2(m,k) = (RMSES3_cs2{k}(m)); 
    end        
    
end


%% PLOTTING 

norows = 3;
nocols = 2;

leadtimes_cs1 = [1:maxleadtime_cs1]*mtu_cs1;
leadtimes_cs2 = [1:maxleadtime_cs2]*mtu_cs2;

leadtimeplot = leadtimes_cs2; 
for g = 1:length(leadtimeplot)
    plotind_cs1(g) = find(leadtimes_cs1 == leadtimeplot(g));
    plotind_cs2(g) = find(leadtimes_cs2 == leadtimeplot(g));
end

xlims = [0 0.64];
ylims = [0 1];
ylimsls = [-0.2 1];
marksz = 6;
lw = 2;

fs=20;

figure
subplot(norows,nocols,1)  %CASE STUDY 1, CRPS 
h1 = errorbar(leadtimes_cs1(plotind_cs1), mean(CRPSS2mean_cs1(:,plotind_cs1)), std(CRPSS2mean_cs1(:,plotind_cs1)), std(CRPSS2mean_cs1(:,plotind_cs1)), 'ko');
hold on 
h2 = errorbar(leadtimes_cs1(plotind_cs1), mean(CRPSS3mean_cs1(:,plotind_cs1)), std(CRPSS3mean_cs1(:,plotind_cs1)), std(CRPSS3mean_cs1(:,plotind_cs1)), 'b^');
grid on 
set(h1, 'MarkerSize', marksz, 'LineWidth', lw)
set(h2, 'MarkerSize', marksz, 'LineWidth', lw)
%xlabel('Forecast Leadtime (MTUs)')
ylabel('CRP Skill Score')
title('Case Study 1')
xlim(xlims)
ylim(ylims)
set(gca, 'YTick', [ylims(1):0.2:ylims(2)])
legend([h1 h2], legnames, 'Orientation', 'Horizontal', 'FontSize', fs)
set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8,'TitleFontSizeMultiplier', 0.9);
set(gca, 'Position', [0.1300    0.7093    0.3347    0.2157])

subplot(norows,nocols,3)  %CASE STUDY 1, LS
h1 = errorbar(leadtimes_cs1(plotind_cs1), mean(LSS2mean_cs1(:,plotind_cs1)), std(LSS2mean_cs1(:,plotind_cs1)), std(LSS2mean_cs1(:,plotind_cs1)), 'ko');
hold on 
h2 = errorbar(leadtimes_cs1(plotind_cs1), mean(LSS3mean_cs1(:,plotind_cs1)), std(LSS3mean_cs1(:,plotind_cs1)), std(LSS3mean_cs1(:,plotind_cs1)), 'b^');
grid on 
set(h1, 'MarkerSize', marksz, 'LineWidth', lw)
set(h2, 'MarkerSize', marksz, 'LineWidth', lw)
%xlabel('Forecast Leadtime (MTUs)')
ylabel('Log Skill Score')
xlim(xlims)
ylim(ylimsls)
set(gca, 'YTick', [ylimsls(1):0.2:ylimsls(2)])
hold on 
plot(xlims, [0 0], 'k')
set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8,'TitleFontSizeMultiplier', 0.9);
set(gca, 'Position', [0.1300    0.4096    0.3347    0.2157])

subplot(norows,nocols,5)  %CASE STUDY 1, RMSES 
h1 = errorbar(leadtimes_cs1(plotind_cs1), mean(RMSES2mean_cs1(:,plotind_cs1)), std(RMSES2mean_cs1(:,plotind_cs1)), std(RMSES2mean_cs1(:,plotind_cs1)), 'ko');
hold on 
h2 = errorbar(leadtimes_cs1(plotind_cs1), mean(RMSES3mean_cs1(:,plotind_cs1)), std(RMSES3mean_cs1(:,plotind_cs1)), std(RMSES3mean_cs1(:,plotind_cs1)), 'b^');
grid on 
set(h1, 'MarkerSize', marksz, 'LineWidth', lw)
set(h2, 'MarkerSize', marksz, 'LineWidth', lw)
%xlabel('Forecast Leadtime (MTUs)')
ylabel('RMSE Skill Score')
xlim(xlims)
ylim(ylims)
set(gca, 'YTick', [ylims(1):0.2:ylims(2)])
set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8)
set(gca, 'Position', [0.1300    0.1100    0.3347    0.2157])

%CASE STUDY 2:
subplot(norows,nocols,2)  %CASE STUDY 2, CRPS 
h1 = errorbar(leadtimes_cs2(plotind_cs2), mean(CRPSS2mean_cs2(:,plotind_cs2)), std(CRPSS2mean_cs2(:,plotind_cs2)), std(CRPSS2mean_cs2(:,plotind_cs2)), 'ko');
hold on 
h2 = errorbar(leadtimes_cs2(plotind_cs2), mean(CRPSS3mean_cs2(:,plotind_cs2)), std(CRPSS3mean_cs2(:,plotind_cs2)), std(CRPSS3mean_cs2(:,plotind_cs2)), 'b^');
grid on 
set(h1, 'MarkerSize', marksz, 'LineWidth', lw)
set(h2, 'MarkerSize', marksz, 'LineWidth', lw)
%xlabel('Forecast Leadtime (MTUs)')
%ylabel('CRP Skill Score')
title('Case Study 2')
xlim(xlims)
ylim(ylims)
set(gca, 'YTick', [ylims(1):0.2:ylims(2)])
set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8,'TitleFontSizeMultiplier', 0.9);
set(gca, 'Position', [0.5703    0.7093    0.3347    0.2157])

subplot(norows,nocols,4)  %CASE STUDY 2, LS
h1 = errorbar(leadtimes_cs2(plotind_cs2), mean(LSS2mean_cs2(:,plotind_cs2)), std(LSS2mean_cs2(:,plotind_cs2)), std(LSS2mean_cs2(:,plotind_cs2)), 'ko');
hold on 
h2 = errorbar(leadtimes_cs2(plotind_cs2), mean(LSS3mean_cs2(:,plotind_cs2)), std(LSS3mean_cs2(:,plotind_cs2)), std(LSS3mean_cs2(:,plotind_cs2)), 'b^');
grid on 
set(h1, 'MarkerSize', marksz, 'LineWidth', lw)
set(h2, 'MarkerSize', marksz, 'LineWidth', lw)
%xlabel('Forecast Leadtime (MTUs)')
%ylabel('Log Skill Score')
xlim(xlims)
ylim(ylims)
set(gca, 'YTick', [ylims(1):0.2:ylims(2)])
hold on 
set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8)
set(gca, 'Position', [0.5703    0.4096    0.3347    0.2157])

subplot(norows,nocols,6)  %CASE STUDY 2, RMSES 
h1 = errorbar(leadtimes_cs2(plotind_cs2), mean(RMSES2mean_cs2(:,plotind_cs2)), std(RMSES2mean_cs2(:,plotind_cs2)), std(RMSES2mean_cs2(:,plotind_cs2)), 'ko');
hold on 
h2 = errorbar(leadtimes_cs2(plotind_cs2), mean(RMSES3mean_cs2(:,plotind_cs2)), std(RMSES3mean_cs2(:,plotind_cs2)), std(RMSES3mean_cs2(:,plotind_cs2)), 'b^');
grid on 
set(h1, 'MarkerSize', marksz, 'LineWidth', lw)
set(h2, 'MarkerSize', marksz, 'LineWidth', lw)
xlabel('Forecast Leadtime (MTUs)')
%ylabel('RMSE Skill Score')
xlim(xlims)
ylim(ylims)
set(gca, 'YTick', [ylims(1):0.2:ylims(2)])
set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8)
set(gca, 'Position', [0.5703    0.1100    0.3347    0.2157])


annotation(gcf,'textbox',...
    [0.129354725787631 0.898177920685959 0.0328389731621937 0.0267952840300107],...
    'String',{'a)'},...
    'LineStyle','none',...
    'FontSize',fs,...
    'FontName','Helvetica',...
    'FitBoxToText','off');

annotation(gcf,'textbox',...
    [0.575095682613767 0.894962486602358 0.0328389731621938 0.0267952840300107],...
    'String','b)',...
    'LineStyle','none',...
    'FontSize',fs,...
    'FontName','Helvetica',...
    'FitBoxToText','off');

annotation(gcf,'textbox',...
    [0.134022170361723 0.594855305466238 0.0328389731621938 0.0267952840300107],...
    'String','c)',...
    'LineStyle','none',...
    'FontSize',fs,...
    'FontName','Helvetica',...
    'FitBoxToText','off');

annotation(gcf,'textbox',...
    [0.575095682613766 0.595927116827438 0.0328389731621938 0.0267952840300107],...
    'String','d)',...
    'LineStyle','none',...
    'FontSize',fs,...
    'FontName','Helvetica',...
    'FitBoxToText','off');

annotation(gcf,'textbox',...
    [0.131688448074675 0.296891747052518 0.0328389731621938 0.0267952840300107],...
    'String','e)',...
    'LineStyle','none',...
    'FontSize',fs,...
    'FontName','Helvetica',...
    'FitBoxToText','off');

annotation(gcf,'textbox',...
    [0.572761960326718 0.296891747052519 0.0328389731621938 0.0267952840300107],...
    'String','f)',...
    'LineStyle','none',...
    'FontSize',fs,...
    'FontName','Helvetica',...
    'FitBoxToText','off');



