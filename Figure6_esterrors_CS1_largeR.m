%Figure: Plot of estimated and true transition errors vs value of state at
%previous time, for both case studies, for a specific obs density.  

clear all
close all

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%USER INPUTS:

%CASE STUDY 1 
%proposed method + truth
load alldata_train_obserr_CS1_PM.mat
alldata_t = alldata_true; alldata_1 = alldata; sampselec_1 = sampselec; noICs_1 = noICs; varsel_1 = varsel; 
f = find(alldata(:,1) <0.04); alldata_1 = alldata_1(f,:);
clearvars -except *_1 *_t decorrtimestep

%method B2 
load alldata_train_obserr_CS1_B2.mat
alldata_3 = alldata; sampselec_3 = sampselec; noICs_3 = noICs; varsel_3 = varsel;
clearvars -except *_1 *_2 *_3 *_t decorrtimestep 

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

%% CASE STUDY 1

norows = 1;
nocols = 3;
pointsize = 10;

fs =  24;

ylims_1 = [-0.08 0.06];
xlims_1 = [-10 15];

figure
%CASE STUDY 1:
subplot(norows,nocols,1)  %the true ones 
plot(alldata_t(:,2), alldata_t(:,1), '.', 'MarkerSize', 10)
xlabel('{\it{\bfx}}_{\itj}_-_1[k]')
ylabel('{\bf{\eta}}_{\itj}[k]')
ylim(ylims_1)
xlim(xlims_1)
title('Truth','FontWeight','Normal')
set(gca, 'YTick', [-0.08:0.02:0.06])
set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8)

subplot(norows,nocols,2)  %estimated ones
plot(alldata_1(:,2), alldata_1(:,1), '.', 'MarkerSize', 10)
ylim(ylims_1)
xlim(xlims_1)
title('Proposed Method','FontWeight','Normal')
set(gca, 'YTick', [-0.08:0.02:0.06])
set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8)

subplot(norows,nocols,3)  %4d-var
plot(alldata_3(:,2), alldata_3(:,1), '.', 'MarkerSize', 10)
ylim(ylims_1)
xlim(xlims_1)
title('B2 Method','FontWeight','Normal')
set(gca, 'YTick', [-0.08:0.02:0.06])
set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8)


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



