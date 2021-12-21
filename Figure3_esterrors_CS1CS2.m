%Figure: Plot of estimated and true transition errors vs value of state at
%previous time, for both case studies, for a specific obs density.  

clear all
close all

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%USER INPUTS:

%CASE STUDY 1 
%proposed method abd true transition errors: 
load alldata_train_CS1_PM.mat
alldata_1 = alldata; alldata_t = alldata_true; sampselec_1 = sampselec; noICs_1 = noICs; varsel_1 = varsel; decorrtimestep1 = decorrtimestep;
clearvars -except *_1 *_t decorrtimestep1 

%method B1
load trainresults_CS1_B1_lambda=1.45.mat
driftall_2 = driftall; xt1all_2 = xt1all; 
clearvars -except *_1 *_2 *_3 *_t decorrtimestep1 

%method B2
load alldata_train_CS1_B2.mat
alldata_3 = alldata; sampselec_3 = sampselec; noICs_3 = noICs; varsel_3 = varsel;
clearvars -except *_1 *_2 *_3 *_t decorrtimestep1 

%CASE STUDY 2 
%proposed method and true transition errors: 
load alldata_train_CS2_PM.mat
alldata_12 = alldata; alldata_t_2 = alldata_true; sampselec_12 = sampselec; noICs_12 = noICs; varsel_12 = varsel; decorrtimestep2 = decorrtimestep;

%method B1 
load trainresults_CS2_B1_lambda=2.mat
driftall_22 = [NaN*ones(size(driftall,1),1), driftall(:,2:end)]; xt1all_22 = [NaN*ones(size(driftall,1),1), xt1all(:,2:end)]; driftall_t1_22 = [NaN*ones(size(driftall,1),1), driftall(:,1:end-1)]; 
clearvars -except *_22 *_12 *_1 *_2 *_3 *_t* decorrtimestep1 decorrtimestep2 sampselec_new 

%method B2 
load alldata_train_CS2_B2.mat
%noICs_32 = noICs; stpt_32 = stpt; %sampselec_32 = sampselec;
%varplot_2 = varsel;
alldata_32 = alldata; sampselec_32 = sampselec; noICs_32 = noICs; varsel_32 = varsel;

clearvars -except *_32 *_22 *_12 *_1 *_2 *_3 *_t* decorrtimestep1 decorrtimestep2 sampselec_new

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

%% CASE STUDY 1

norows = 2;
nocols = 4;
pointsize = 10;

fs = 20;  %fontsize 

ylims_1 = [-0.08 0.08];
xlims_1 = [-10 15];
yticks_1 = [-0.08:0.04:0.08];
xticks = [-10:5:15];

%first extract the same points in the analysis increment one:
lengtheach_2 = reshape([1:size(driftall_2,2)], size(driftall_2,2)/noICs_1, noICs_1)';
keepinds_2 = lengtheach_2(:, sampselec_1-1);   %-1 index needed because first value ignored in the ETKF analysis increment data. 
keepinds_2 = keepinds_2(:);

figure
%CASE STUDY 1:
subplot(norows,nocols,1)  %the true ones 
plot(alldata_t(:,2), alldata_t(:,1), '.')
ylabel('{\bf{\eta}}_{\itj}[k]')
ylim(ylims_1)
xlim(xlims_1)
title('Truth','FontWeight','Normal')
set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8,'YTick', yticks_1, 'XTick', xticks)

subplot(norows,nocols,2)  %estimated ones
plot(alldata_1(:,2), alldata_1(:,1), '.')
ylim(ylims_1)
xlim(xlims_1)
title('Proposed Method','FontWeight','Normal')
set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8,'YTick', yticks_1, 'XTick', xticks)

subplot(norows,nocols,3)  %4d-var
plot(alldata_3(:,2), alldata_3(:,1), '.')
ylim(ylims_1)
xlim(xlims_1)
title('B2 Method','FontWeight','Normal')
set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8,'YTick', yticks_1, 'XTick', xticks)
xlabel('{\it{\bfx}}_{\itj}_-_1[k]')

subplot(norows,nocols,4)  %from analysis increments 
plot(xt1all_2(varsel_1, keepinds_2), driftall_2(varsel_1,keepinds_2), '.')
ylim(ylims_1)
xlim(xlims_1)
title('B1 Method','FontWeight','Normal')
set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8,'YTick', yticks_1, 'XTick', xticks)

%annotation:
annotation(gcf,'textbox',...
    [0.132836499712148 0.861742424242425 0.0404507772020725 0.0549242424242432],...
    'String','a)',...
    'LineStyle','none',...
    'FontSize',fs,...
    'FontName','Helvetica',...
    'FitBoxToText','off');

annotation(gcf,'textbox',...
    [0.338363270005759 0.863636363636366 0.0657818077144484 0.0549242424242432],...
    'String','b)',...
    'LineStyle','none',...
    'FontSize',fs,...
    'FontName','Helvetica',...
    'FitBoxToText','off');

annotation(gcf,'textbox',...
    [0.543890040299371 0.861742424242428 0.0657818077144484 0.0549242424242432],...
    'String','c)',...
    'LineStyle','none',...
    'FontSize',fs,...
    'FontName','Helvetica',...
    'FitBoxToText','off');

annotation(gcf,'textbox',...
    [0.751143926309739 0.865530303030307 0.0657818077144484 0.0549242424242432],...
    'String','d)',...
    'LineStyle','none',...
    'FontSize',fs,...
    'FontName','Helvetica',...
    'FitBoxToText','off');


%% CASE STUDY 2 

ylims_2 = [-10 15];
xlims_2 = [-10 15];
yticks = [-10:5:15];
xticks = [-10:5:15];

%first extract the same points in the analysis increment one:
lengtheach_22 = reshape([1:size(driftall_22,2)], size(driftall_22,2)/noICs_12, noICs_12)';
keepinds_22 = lengtheach_22(:, sampselec_12);
keepinds_22 = keepinds_22(:);

%first create the interpolating colours:
x12 = alldata_12(:,2)'; %tx_estall_12(varsel_12,:);   %proposed method 
y12 = alldata_12(:,3)';%tx_estall_12(varsel_12-1,:);
z12 = alldata_12(:,1)';%errx_estall_12(varsel_12,:);

x32 = alldata_32(:,2)'; %tx_estall_32(varsel_12,:);   %method B2
y32 = alldata_32(:,3)'; %tx_estall_32(varsel_12-1,:);
z32 = alldata_32(:,1)'; %errx_estall_32(varsel_12,:);

xt2 = alldata_t_2(:,2)'; %tx_t_estall_12(varsel_12,:);   %truth 
yt2 = alldata_t_2(:,3)'; %tx_t_estall_12(varsel_12-1,:);
zt2 = alldata_t_2(:,1)'; %errx_t_estall_12(varsel_12,:);

x22 = xt1all_22(varsel_12,keepinds_22);   %method B1
y22 = xt1all_22(varsel_12-1,keepinds_22);
z22 = driftall_22(varsel_12,keepinds_22);

figure
scatter([xt2, x12, x32], [yt2, y12, y32], pointsize, [zt2, z12, z32]);
cm = colormap('jet');
zm = linspace(min([zt2,z12,z32])-0.005, max([zt2,z12,z32])+0.005, size(cm,1)+1); 
close 

%now apply into each bin:
for l = 1:length(zt2)
    tmp = find(zm > zt2(l));
    Ct(l) = tmp(1)-1;  %color bin 
end

for l = 1:length(z12)
    tmp = find(zm > z12(l));
    C12(l) = tmp(1)-1;  %color bin 
end

for l = 1:length(z22)
    tmp = find(zm > z22(l));
    if isempty(tmp)
        %larger - just make equal to maximum:
        C22(l) = length(zm)-1;
    elseif tmp(1) == 1
        C22(l) = 1;  %just make equal to min 
    else        
        C22(l) = tmp(1)-1;  %color bin 
    end
end

for l = 1:length(z32)
    tmp = find(zm > z32(l));
    C32(l) = tmp(1)-1;  %color bin 
end

colticks = round(linspace(zm(2), zm(end), 11), 2);

%first define rectangle bounds for annotation
rec1_y = [5 12];  rec1_x = [-6 2];
rec2_y = [-1.1 4]; rec2_x = [1.5 7.5];

rectx = [rec1_x; rec2_x];  recty = [rec1_y; rec2_y];
rectx = rectx(1,:); recty = recty(1,:); 

subplot(norows,nocols,5)  %the true ones 
scatter(xt2, yt2, pointsize, cm(Ct,:));
ylabel('{\it{\bfx}}_{\itj}_-_1[k-1]')
ylim(ylims_2)
xlim(xlims_2)
set(gca, 'YTick', yticks, 'XTick', xticks)
set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8)
box on
colorbar('Ticks',[0:0.1:1],'TickLabels',colticks, 'Location', 'southoutside')

annotation(gcf,'textbox',...
    [0.506208333333333 0.439834024896266 0.0406666666666665 0.0373443983402489],...
    'String',{'{\bf{\eta}}_{\itj}[k]'},...
    'LineStyle','none',...
    'FontSize',fs,...
    'FontName','Helvetica',...
    'FitBoxToText','off');


subplot(norows,nocols,6)  %proposed method 
scatter(x12, y12, pointsize, cm(C12,:));
ylim(ylims_2)
xlim(xlims_2)
set(gca, 'YTick', yticks, 'XTick', xticks)
set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8)
box on


subplot(norows,nocols,7)  %method B2
scatter(x32, y32, pointsize, cm(C32,:));
ylim(ylims_2)
xlim(xlims_2)
set(gca, 'YTick', yticks, 'XTick', xticks)
xlabel('{\it{\bfx}}_{\itj}_-_1[k]')
set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8)
box on 


subplot(norows,nocols,8)  %method B1 
scatter(x22, y22, pointsize, cm(C22,:));
ylim(ylims_2)
xlim(xlims_2)
set(gca, 'YTick', yticks, 'XTick', xticks)
set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8)
box on 

%now reposition axes:
subplot(norows, nocols,5)
set(gca, 'Position',[0.13 0.164979253112033 0.156648936170213 0.341162790697674]);
subplot(norows, nocols,6)
set(gca, 'Position',[0.336117021276596 0.164979253112033 0.156648936170213 0.341162790697674]);
subplot(norows, nocols,7)
set(gca, 'Position',[0.542234042553191 0.164979253112033 0.156648936170213 0.341162790697674]);
subplot(norows, nocols,8)
set(gca, 'Position',[0.748351063829787 0.164979253112033 0.156648936170213 0.341162790697674]);


annotation(gcf,'textbox',...
    [0.131794833045481 0.441617942914624 0.0404507772020725 0.0549242424242432],...
    'String','e)',...
    'LineStyle','none',...
    'FontSize',fs,...
    'FontName','Helvetica',...
    'FitBoxToText','off');

annotation(gcf,'textbox',...
    [0.337523999712148 0.440580598516284 0.0404507772020725 0.0549242424242432],...
    'String','f)',...
    'LineStyle','none',...
    'FontSize',fs,...
    'FontName','Helvetica',...
    'FitBoxToText','off');

annotation(gcf,'textbox',...
    [0.544294833045481 0.440580598516284 0.0404507772020724 0.0549242424242432],...
    'String','g)',...
    'LineStyle','none',...
    'FontSize',fs,...
    'FontName','Helvetica',...
    'FitBoxToText','off');

annotation(gcf,'textbox',...
    [0.749503166378814 0.439543254117944 0.0404507772020725 0.0549242424242432],...
    'String','h)',...
    'LineStyle','none',...
    'FontSize',fs,...
    'FontName','Helvetica',...
    'FitBoxToText','off');


