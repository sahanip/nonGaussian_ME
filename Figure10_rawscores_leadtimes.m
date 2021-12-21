%Figure: raw scores for different lead times for all 3 methods
%method for Case study 1.  Plots CRPS, LS, RMSE.

clear all
close all

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%USER INPUTS:

%CASE STUDY 1:
%following is to load obsind and hidind:
load assimilateresults_CS1_PM.mat
obsind_cs1 = obsind; hidind_cs1 = hidind; x_true_cs1 = x_true;
clearvars -except *_cs1 

%load proposed approach 
load forecastresults_CS1_PM.mat
CRPSall1_cs1 = CRPSall; ENSVARall1_cs1 = ENSVARall; LSall1_cs1 = LSall; ENSMEANERRall1_cs1 = ENSMEANERRall;
clearvars -except *_cs1

%now method b1 
load forecastresults_CS1_B1_alpha=0.8.mat
CRPSall2_cs1 = CRPSall; ENSVARall2_cs1 = ENSVARall; LSall2_cs1 = LSall; ENSMEANERRall2_cs1 = ENSMEANERRall;
clearvars -except *_cs1

%Now method b2 
%load my results scores:
load forecastresults_CS1_B2.mat
CRPSall3_cs1 = CRPSall; ENSVARall3_cs1 = ENSVARall; LSall3_cs1 = LSall; ENSMEANERRall3_cs1 = ENSMEANERRall;
maxleadtime_cs1 = maxleadtime;
clearvars -except *_cs1 


%CASE STUDY 2:
%load proposed approach 
load forecastresults_CS2_PM.mat
CRPSall1_cs2 = CRPSall; ENSVARall1_cs2 = ENSVARall; LSall1_cs2 = LSall; ENSMEANERRall1_cs2 = ENSMEANERRall;
clearvars -except *_cs1 *_cs2 

%now method b1 
load forecastresults_CS2_B1_alpha=0.8.mat
CRPSall2_cs2 = CRPSall; ENSVARall2_cs2 = ENSVARall; LSall2_cs2 = LSall; ENSMEANERRall2_cs2 = ENSMEANERRall;
clearvars -except *_cs1 *_cs2

%Now method b2 
load forecastresults_CS2_B2.mat
CRPSall3_cs2 = CRPSall; ENSVARall3_cs2 = ENSVARall; LSall3_cs2 = LSall; ENSMEANERRall3_cs2 = ENSMEANERRall;
maxleadtime_cs2 = maxleadtime;
clearvars -except *_cs1 *_cs2 

mtu_cs1 = 0.02;  %mtu of spacing of observations 
mtu_cs2 = 0.04;  %mtu of spacing of observations

legnames = {'Proposed Method', 'B2 Method', 'B1 Method', }; %order: 1,3,2

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

noICs = length(CRPSall1_cs1);

%% CASE STUDY 1:

%calculate the mean score across observed and hidden variables, across
%lead times:

CRPS1_LT_cs1 = NaN*ones(1, maxleadtime_cs1);   
LS1_LT_cs1 = CRPS1_LT_cs1; RMSE1_LT_cs1 = CRPS1_LT_cs1;
CRPS2_LT_cs1 = CRPS1_LT_cs1;  LS2_LT_cs1 = CRPS1_LT_cs1;  RMSE2_LT_cs1 = CRPS1_LT_cs1; 
CRPS3_LT_cs1 = CRPS1_LT_cs1;  LS3_LT_cs1 = CRPS1_LT_cs1;  RMSE3_LT_cs1 = CRPS1_LT_cs1; 

for m = 1:maxleadtime_cs1 
    crpstemp1_all_cs1 = []; crpstemp2_all_cs1 = []; crpstemp3_all_cs1 = [];  lstemp1_all_cs1 = []; lstemp2_all_cs1 = [];  lstemp3_all_cs1 = [];  rmsetemp1_all_cs1 = []; rmsetemp2_all_cs1 = [];  rmsetemp3_all_cs1 = [];
    
    for j = 1:noICs
        
        %crps:
        crpstemp1_cs1 = CRPSall1_cs1{j};
        crpstemp1_all_cs1 = [crpstemp1_all_cs1, crpstemp1_cs1(:,:,m)];
        
        crpstemp2_cs1 = CRPSall2_cs1{j};
        crpstemp2_all_cs1 = [crpstemp2_all_cs1, crpstemp2_cs1(:,:,m)];
        
        crpstemp3_cs1 = CRPSall3_cs1{j};
        crpstemp3_all_cs1 = [crpstemp3_all_cs1, crpstemp3_cs1(:,:,m)];
        
        %log score:
        lstemp1_cs1 = LSall1_cs1{j};
        lstemp1_all_cs1 = [lstemp1_all_cs1, lstemp1_cs1(:,:,m)];
        
        lstemp2_cs1 = LSall2_cs1{j};
        lstemp2_all_cs1 = [lstemp2_all_cs1, lstemp2_cs1(:,:,m)];
        
        lstemp3_cs1 = LSall3_cs1{j};
        lstemp3_all_cs1 = [lstemp3_all_cs1, lstemp3_cs1(:,:,m)];
        
        %rmse:
        %have to collate error in ensemble mean and then calculate RMS at
        %end.
        rmsetemp1_cs1 = sqrt(nanmean((ENSMEANERRall1_cs1{j}).^2,2));  %just calculate rmse 
        rmsetemp1_all_cs1 = [rmsetemp1_all_cs1, rmsetemp1_cs1(:,:,m)]; 
        
        rmsetemp2_cs1 = sqrt(nanmean((ENSMEANERRall2_cs1{j}).^2,2));  %just calculate rmse 
        rmsetemp2_all_cs1 = [rmsetemp2_all_cs1, rmsetemp2_cs1(:,:,m)]; 
        
        rmsetemp3_cs1 = sqrt(nanmean((ENSMEANERRall3_cs1{j}).^2,2));  %just calculate rmse 
        rmsetemp3_all_cs1 = [rmsetemp3_all_cs1, rmsetemp3_cs1(:,:,m)];
    end
    
    %crps:
    CRPS1_LT_cs1(1,m) = nanmean(crpstemp1_all_cs1(:)); 
    CRPS3_LT_cs1(1,m) = nanmean(crpstemp3_all_cs1(:)); 
    CRPS2_LT_cs1(1,m) = nanmean(crpstemp2_all_cs1(:)); 
    
    %ls:
    LS1_LT_cs1(1,m) = nanmean(lstemp1_all_cs1(:)); 
    LS3_LT_cs1(1,m) = nanmean(lstemp3_all_cs1(:)); 
    LS2_LT_cs1(1,m) = nanmean(lstemp2_all_cs1(:)); 
    
    %rmse:     
    RMSE1_LT_cs1(1,m) = nanmean(rmsetemp1_all_cs1(:)); 
    RMSE3_LT_cs1(1,m) = nanmean(rmsetemp3_all_cs1(:)); 
    RMSE2_LT_cs1(1,m) = nanmean(rmsetemp2_all_cs1(:)); 
end

%% NOW FOR CASE STUDY 2:

%calculate the mean score across observed and hidden variables, across
%lead times:

CRPS1_LT_cs2 = NaN*ones(1, maxleadtime_cs2);   
LS1_LT_cs2 = CRPS1_LT_cs2; RMSE1_LT_cs2 = CRPS1_LT_cs2;
CRPS2_LT_cs2 = CRPS1_LT_cs2;  LS2_LT_cs2 = CRPS1_LT_cs2;  RMSE2_LT_cs2 = CRPS1_LT_cs2; 
CRPS3_LT_cs2 = CRPS1_LT_cs2;  LS3_LT_cs2 = CRPS1_LT_cs2;  RMSE3_LT_cs2 = CRPS1_LT_cs2; 

for m = 1:maxleadtime_cs2 
    crpstemp1_all_cs2 = []; crpstemp2_all_cs2 = []; crpstemp3_all_cs2 = [];  lstemp1_all_cs2 = []; lstemp2_all_cs2 = [];  lstemp3_all_cs2 = [];  rmsetemp1_all_cs2 = []; rmsetemp2_all_cs2 = [];  rmsetemp3_all_cs2 = [];
    
    for j = 1:noICs
        
        %crps:
        crpstemp1_cs2 = CRPSall1_cs2{j};
        crpstemp1_all_cs2 = [crpstemp1_all_cs2, crpstemp1_cs2(:,:,m)];
        
        crpstemp2_cs2 = CRPSall2_cs2{j};
        crpstemp2_all_cs2 = [crpstemp2_all_cs2, crpstemp2_cs2(:,:,m)];
        
        crpstemp3_cs2 = CRPSall3_cs2{j};
        crpstemp3_all_cs2 = [crpstemp3_all_cs2, crpstemp3_cs2(:,:,m)];
        
        %log score:
        lstemp1_cs2 = LSall1_cs2{j};
        lstemp1_all_cs2 = [lstemp1_all_cs2, lstemp1_cs2(:,:,m)];
        
        lstemp2_cs2 = LSall2_cs2{j};
        lstemp2_all_cs2 = [lstemp2_all_cs2, lstemp2_cs2(:,:,m)];
        
        lstemp3_cs2 = LSall3_cs2{j};
        lstemp3_all_cs2 = [lstemp3_all_cs2, lstemp3_cs2(:,:,m)];
        
        %rmse:
        %have to collate error in ensemble mean and then calculate RMS at
        %end.
        rmsetemp1_cs2 = sqrt(nanmean((ENSMEANERRall1_cs2{j}).^2,2));  %just calculate rmse 
        rmsetemp1_all_cs2 = [rmsetemp1_all_cs2, rmsetemp1_cs2(:,:,m)]; 
        
        rmsetemp2_cs2 = sqrt(nanmean((ENSMEANERRall2_cs2{j}).^2,2));  %just calculate rmse 
        rmsetemp2_all_cs2 = [rmsetemp2_all_cs2, rmsetemp2_cs2(:,:,m)]; 
        
        rmsetemp3_cs2 = sqrt(nanmean((ENSMEANERRall3_cs2{j}).^2,2));  %just calculate rmse 
        rmsetemp3_all_cs2 = [rmsetemp3_all_cs2, rmsetemp3_cs2(:,:,m)];
    end
    
    %crps:
    CRPS1_LT_cs2(1,m) = nanmean(crpstemp1_all_cs2(:)); 
    CRPS3_LT_cs2(1,m) = nanmean(crpstemp3_all_cs2(:)); 
    CRPS2_LT_cs2(1,m) = nanmean(crpstemp2_all_cs2(:)); 
    
    %ls:
    LS1_LT_cs2(1,m) = nanmean(lstemp1_all_cs2(:)); 
    LS3_LT_cs2(1,m) = nanmean(lstemp3_all_cs2(:)); 
    LS2_LT_cs2(1,m) = nanmean(lstemp2_all_cs2(:)); 
    
    %rmse:     
    RMSE1_LT_cs2(1,m) = nanmean(rmsetemp1_all_cs2(:)); 
    RMSE3_LT_cs2(1,m) = nanmean(rmsetemp3_all_cs2(:)); 
    RMSE2_LT_cs2(1,m) = nanmean(rmsetemp2_all_cs2(:)); 
end

%% PLOTTING

norows = 3;
nocols = 2;

leadtimes_cs1 = [1:maxleadtime_cs1]*mtu_cs1; 
leadtimes_cs2 = [1:maxleadtime_cs2]*mtu_cs2; 

lw = 2  %linewidth 
ms = 16   %markersize 

xtickpts = [0:0.1:0.6];

fs=20 

%first for CS1:
figure
subplot(norows, nocols, 1)  %CRPS
h21 = plot(leadtimes_cs1, CRPS1_LT_cs1, '.-r', leadtimes_cs1, CRPS3_LT_cs1, '.-b',leadtimes_cs1, CRPS2_LT_cs1, '.-k')
set(h21, 'LineWidth', lw, 'MarkerSize', ms);
legend(legnames, 'Orientation', 'Horizontal', 'FontSize', fs)
ylabel('CRPS')
%set(gca, 'XTick', leadtimes_cs1(1:2:end))
set(gca, 'XTick', xtickpts)
set(gca,'XMinorTick','on')
xlim([0 (1+ maxleadtime_cs1)*mtu_cs1])
grid on 
%xlabel('Forecast Leadtime (MTUs)')
title('Case Study 1')
set(gca, 'Position', [0.1300    0.7093    0.3347    0.2157],'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8,'TitleFontSizeMultiplier', 0.9);


subplot(norows,nocols,3)
h21 = plot(leadtimes_cs1, LS1_LT_cs1, '.-r', leadtimes_cs1, LS3_LT_cs1, '.-b', leadtimes_cs1, LS2_LT_cs1, '.-k');
set(h21, 'LineWidth', lw, 'MarkerSize', ms);
ylabel('LS')
set(gca, 'XTick', leadtimes_cs1)
grid on 
xlim([0 (1+ maxleadtime_cs1)*mtu_cs1])
set(gca, 'XTick', xtickpts)
set(gca,'XMinorTick','on')
set(gca, 'Position', [0.1300    0.4096    0.3347    0.2157],'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8)


subplot(norows,nocols,5)
h21 = plot(leadtimes_cs1, RMSE1_LT_cs1, '.-r', leadtimes_cs1, RMSE3_LT_cs1, '.-b', leadtimes_cs1, RMSE2_LT_cs1, '.-k');
set(h21, 'LineWidth', lw, 'MarkerSize', ms);
set(gca, 'XTick', leadtimes_cs1)
grid on 
ylabel('RMSE')
xlim([0 (1+ maxleadtime_cs1)*mtu_cs1])
set(gca, 'XTick', xtickpts)
set(gca,'XMinorTick','on')
set(gca, 'Position', [0.1300    0.1100    0.3347    0.2157],'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8)

%CASE STUDY 2:
subplot(norows, nocols, 2)  %CRPS
h21 = plot(leadtimes_cs2, CRPS1_LT_cs2, '.-r', leadtimes_cs2, CRPS3_LT_cs2, '.-b',leadtimes_cs2, CRPS2_LT_cs2, '.-k')
set(h21, 'LineWidth', lw, 'MarkerSize', ms);
set(gca, 'XTick', xtickpts)
set(gca,'XMinorTick','on')
xlim([0 (1+ maxleadtime_cs2)*mtu_cs2])
grid on 
title('Case Study 2')
set(gca, 'Position', [0.5703    0.7093    0.3347    0.2157],'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8,'TitleFontSizeMultiplier', 0.9);

subplot(norows,nocols,4)
h21 = plot(leadtimes_cs2, LS1_LT_cs2, '.-r', leadtimes_cs2, LS3_LT_cs2, '.-b', leadtimes_cs2, LS2_LT_cs2, '.-k');
set(h21, 'LineWidth', lw, 'MarkerSize', ms);
set(gca, 'XTick', leadtimes_cs2)
grid on 
xlim([0 (1+ maxleadtime_cs2)*mtu_cs2])
set(gca, 'XTick', xtickpts)
set(gca,'XMinorTick','on')
set(gca, 'Position', [0.5703    0.4096    0.3347    0.2157],'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8)

subplot(norows,nocols,6)
h21 = plot(leadtimes_cs2, RMSE1_LT_cs2, '.-r', leadtimes_cs2, RMSE3_LT_cs2, '.-b', leadtimes_cs2, RMSE2_LT_cs2, '.-k');
set(h21, 'LineWidth', lw, 'MarkerSize', ms);
set(gca, 'XTick', leadtimes_cs2)
grid on 
xlim([0 (1+ maxleadtime_cs2)*mtu_cs2])
set(gca, 'XTick', xtickpts)
set(gca,'XMinorTick','on')
xlabel('Forecast Leadtime (MTUs)')
set(gca, 'Position', [0.5703    0.1100    0.3347    0.2157],'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8)


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








