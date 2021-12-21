%Figure: ensemble spread vs rmse plot

clear all
close all

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%USER INPUTS:

%CASE STUDY 1:
%proposed method 
load forecastresults_CS1_PM.mat
ENSVARall1_cs1 = ENSVARall; ENSMEANERRall1_cs1 = ENSMEANERRall;
clearvars -except *_cs1 

%Now method b1 
load forecastresults_CS1_B1_alpha=0.8.mat
ENSVARall2_cs1 = ENSVARall; ENSMEANERRall2_cs1 = ENSMEANERRall;
clearvars -except  *_cs1 

%method b2 
load forecastresults_CS1_B2.mat
ENSVARall3_cs1 = ENSVARall; ENSMEANERRall3_cs1 = ENSMEANERRall;
clearvars -except *_cs1 

%-----------------
%CASE STUDY 2:
%proposed method 
load forecastresults_CS2_PM.mat 
ENSVARall1_cs2 = ENSVARall; ENSMEANERRall1_cs2 = ENSMEANERRall;
clearvars -except *_cs2 *_cs1

%Now method b1
load forecastresults_CS2_B1_alpha=0.8.mat
ENSVARall2_cs2 = ENSVARall; ENSMEANERRall2_cs2 = ENSMEANERRall;
clearvars -except  *_cs2 *_cs1

%Now method b1 
load forecastresults_CS2_B2.mat 
ENSVARall3_cs2 = ENSVARall; ENSMEANERRall3_cs2 = ENSMEANERRall;
clearvars -except *_cs2 *_cs1

nobins = 10;  %for binning the ensemble variance

legnames = {'Proposed Approach', 'B2 Method','B1 Method'};  %order: 1, 3, 2

leadtimecons_cs1 = 30; %FOR CASE STUDY 1
leadtimecons_cs2 = 15; %FOR CASE STUDY 2


%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


%%  CASE STUDY 1:

noICs = length(ENSVARall1_cs1);

ev1_cs1 = []; ev2_cs1 = []; ev3_cs1 = []; 
me1_cs1 = []; me2_cs1 = []; me3_cs1 = []; 

for m = 1:noICs
 
    ensvar1 = ENSVARall1_cs1{m}(:,:,leadtimecons_cs1); ensvar2 = ENSVARall2_cs1{m}(:,:,leadtimecons_cs1); ensvar3 = ENSVARall3_cs1{m}(:,:,leadtimecons_cs1);
    
    %timesel = [2:size(ensvar1,2)];
    timesel = find(~isnan(ensvar1(1,:)));
    
    %First separate into hidden and obs vars:
    ev1_cs1= [ev1_cs1, ensvar1(:,timesel)];  
    ev2_cs1 = [ev2_cs1, ensvar2(:,timesel)];  
    ev3_cs1 = [ev3_cs1, ensvar3(:,timesel)];  

    ensmeanerr1 = (ENSMEANERRall1_cs1{m}(:,:,leadtimecons_cs1)).^2; ensmeanerr2 = (ENSMEANERRall2_cs1{m}(:,:,leadtimecons_cs1)).^2; ensmeanerr3 = (ENSMEANERRall3_cs1{m}(:,:,leadtimecons_cs1)).^2;  %HAVE TO SQUARE BECAUSE CONTAINS JUST MEAN DIFFERENCE! 
    
    me1_cs1 = [me1_cs1, ensmeanerr1(:,timesel)];  
    me2_cs1 = [me2_cs1, ensmeanerr2(:,timesel)];  
    me3_cs1 = [me3_cs1, ensmeanerr3(:,timesel)];  
 
end

ev1_cs1 = ev1_cs1(:); ev2_cs1 = ev2_cs1(:); ev3_cs1 = ev3_cs1(:);
me1_cs1 = me1_cs1(:); me2_cs1 = me2_cs1(:); me3_cs1 = me3_cs1(:);

%now bin ensemble variance into nbins equally populated bins:
[~, a1_cs1] = sort(ev1_cs1); 
[~, a2_cs1] = sort(ev2_cs1);  
[~, a3_cs1] = sort(ev3_cs1);  

binind_cs1 = reshape([1:length(a1_cs1)], length(a1_cs1)/nobins, nobins);  %consistent for all 

%now calculate mean in each bin:
for k = 1:size(binind_cs1,2)
    errbin1_cs1(k) = sqrt(mean(me1_cs1(a1_cs1(binind_cs1(:,k)))));  %rmse in given bin 
    spreadbin1_cs1(k) = sqrt(mean(ev1_cs1(a1_cs1(binind_cs1(:,k)))));   %spread estimate in given bin 
    
    errbin2_cs1(k) = sqrt(mean(me2_cs1(a2_cs1(binind_cs1(:,k)))));  %rmse in given bin 
    spreadbin2_cs1(k) = sqrt(mean(ev2_cs1(a2_cs1(binind_cs1(:,k)))));   %spread estimate in given bin 
    
    errbin3_cs1(k) = sqrt(mean(me3_cs1(a3_cs1(binind_cs1(:,k)))));  %rmse in given bin 
    spreadbin3_cs1(k) = sqrt(mean(ev3_cs1(a3_cs1(binind_cs1(:,k)))));   %spread estimate in given bin 
end


%% CASE STUDY 2:
noICs = length(ENSVARall1_cs2);

ev1_cs2 = []; ev2_cs2 = []; ev3_cs2 = []; 
me1_cs2 = []; me2_cs2 = []; me3_cs2 = []; 

for m = 1:noICs
 
    ensvar1 = ENSVARall1_cs2{m}(:,:,leadtimecons_cs2); ensvar2 = ENSVARall2_cs2{m}(:,:,leadtimecons_cs2); ensvar3 = ENSVARall3_cs2{m}(:,:,leadtimecons_cs2);
    
    %timesel = [2:size(ensvar1,2)];
    timesel = find(~isnan(ensvar1(1,:)));
    
    %First separate into hidden and obs vars:
    ev1_cs2= [ev1_cs2, ensvar1(:,timesel)];  
    ev2_cs2 = [ev2_cs2, ensvar2(:,timesel)];  
    ev3_cs2 = [ev3_cs2, ensvar3(:,timesel)];  

    ensmeanerr1 = (ENSMEANERRall1_cs2{m}(:,:,leadtimecons_cs2)).^2; ensmeanerr2 = (ENSMEANERRall2_cs2{m}(:,:,leadtimecons_cs2)).^2; ensmeanerr3 = (ENSMEANERRall3_cs2{m}(:,:,leadtimecons_cs2)).^2;  %HAVE TO SQUARE BECAUSE CONTAINS JUST MEAN DIFFERENCE! 
    
    me1_cs2 = [me1_cs2, ensmeanerr1(:,timesel)];  
    me2_cs2 = [me2_cs2, ensmeanerr2(:,timesel)];  
    me3_cs2 = [me3_cs2, ensmeanerr3(:,timesel)];  
 
end

ev1_cs2 = ev1_cs2(:); ev2_cs2 = ev2_cs2(:); ev3_cs2 = ev3_cs2(:);
me1_cs2 = me1_cs2(:); me2_cs2 = me2_cs2(:); me3_cs2 = me3_cs2(:);

%now bin ensemble variance into nbins equally populated bins:
[~, a1_cs2] = sort(ev1_cs2); 
[~, a2_cs2] = sort(ev2_cs2);  
[~, a3_cs2] = sort(ev3_cs2);  

binind_cs2 = reshape([1:length(a1_cs2)], length(a1_cs2)/nobins, nobins);  %consistent for all 

%now calculate mean in each bin:
for k = 1:size(binind_cs2,2)
    errbin1_cs2(k) = sqrt(mean(me1_cs2(a1_cs2(binind_cs2(:,k)))));  %rmse in given bin 
    spreadbin1_cs2(k) = sqrt(mean(ev1_cs2(a1_cs2(binind_cs2(:,k)))));   %spread estimate in given bin 
    
    errbin2_cs2(k) = sqrt(mean(me2_cs2(a2_cs2(binind_cs2(:,k)))));  %rmse in given bin 
    spreadbin2_cs2(k) = sqrt(mean(ev2_cs2(a2_cs2(binind_cs2(:,k)))));   %spread estimate in given bin 
    
    errbin3_cs2(k) = sqrt(mean(me3_cs2(a3_cs2(binind_cs2(:,k)))));  %rmse in given bin 
    spreadbin3_cs2(k) = sqrt(mean(ev3_cs2(a3_cs2(binind_cs2(:,k)))));   %spread estimate in given bin 
end

%% PLOT 

norows = 1;
nocols = 2;

fs=24;

lw = 2;
ms = 6;


figure
subplot(norows, nocols,1)  %CASE STUDY 1
h1 = plot(spreadbin1_cs1, errbin1_cs1, 'ro', spreadbin3_cs1, errbin3_cs1, 'b^', spreadbin2_cs1, errbin2_cs1, 'ks'); %NOTE ORDER TO CORRESPOND WITH LEGEND!
set(h1(1), 'MarkerFaceColor','r','MarkerSize', ms)
set(h1(2), 'MarkerFaceColor','b','MarkerSize', ms)
set(h1(3), 'MarkerFaceColor','k','MarkerSize', ms)
legend(legnames, 'Orientation', 'Horizontal', 'Location', 'southoutside', 'FontSize', fs)
ylims = ylim; xlims = xlim;
hold on 
h2 = plot([0, max(ylims(2), xlims(2))], [0, max(ylims(2), xlims(2))], '--k');
set(h2, 'LineWidth', lw)
xlabel('r.m.s. spread'); ylabel('r.m.s. error')
set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8,'TitleFontSizeMultiplier', 0.9);
grid on 
title('Case Study 1')
set(gca, 'Position', [0.1300    0.2211    0.3347    0.7030])

subplot(norows, nocols,2)  %CASE STUDY 2
h1 = plot(spreadbin1_cs2, errbin1_cs2, 'ro', spreadbin3_cs2, errbin3_cs2, 'b^', spreadbin2_cs2, errbin2_cs2, 'ks'); %NOTE ORDER TO CORRESPOND WITH LEGEND!
set(h1(1), 'MarkerFaceColor','r','MarkerSize', ms)
set(h1(2), 'MarkerFaceColor','b','MarkerSize', ms)
set(h1(3), 'MarkerFaceColor','k','MarkerSize', ms)
ylims = ylim; xlims = xlim;
hold on 
h2 = plot([0, max(ylims(2), xlims(2))], [0, max(ylims(2), xlims(2))], '--k');
set(h2, 'LineWidth', lw)
%sxlabel('r.m.s. spread'); ylabel('r.m.s. error')
set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8,'TitleFontSizeMultiplier', 0.9);
grid on 
title('Case Study 2')
set(gca, 'Position', [0.5703    0.2211    0.3347    0.7030])

annotation(gcf,'textbox',...
    [0.13270272812794 0.851318944844124 0.0347478833490122 0.0647482014388473],...
    'String',{'a)'},...
    'LineStyle','none',...
    'FontSize',fs,...
    'FontName','Helvetica',...
    'FitBoxToText','off');

annotation(gcf,'textbox',...
    [0.572025399811854 0.853717026378894 0.0347478833490121 0.0647482014388473],...
    'String','b)',...
    'LineStyle','none',...
    'FontSize',fs,...
    'FontName','Helvetica',...
    'FitBoxToText','off');






