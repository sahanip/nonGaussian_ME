%Figure: this figure shows the benefit of accounting for the non-gaussian
%transition densities for an observed variable

clear all
close all

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%USER INPUTS:
%CASE STUDY 1:
%proposed approach 
load assimilateresults_CS1_PM.mat
ensstates1_1 = ensstates; statesA1_1 = statesA;  tempensstates1_1 = tempensstates; eysample1_1 = eysample;
clearvars -except *_1  

%method B1 
load assimilateresults_CS1_B1_alpha=0.8.mat
ensstates2_1 = ensstates; statesA2_1 = statesA;  tempensstates2_1 = tempensstates; 
clearvars -except *_1  

%method B2 
load assimilateresults_CS1_B2.mat 
ensstates3_1 = ensstates; statesA3_1 = statesA;  tempensstates3_1 = tempensstates; eysample3_1 = eysample; x_true_1 = x_true;
clearvars -except *_1  

legnames = {'Proposed Method', 'B2 Method','B1 Method'};  

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

%% PLOT 

norows = 1;
nocols = 3;

fs=24;
%CS1:
varplot_1 = 9;
timeplot_1 = 82;
adderr2_1 = ensstates2_1 - tempensstates2_1;  %added error for benchmark

grycol = [0.501960784313725 0.501960784313725 0.501960784313725];

%fit probability density:
[y1, x1] = ksdensity(tempensstates1_1(varplot_1,:,timeplot_1));
[y2, x2] = ksdensity(eysample1_1(varplot_1,:,timeplot_1-1));  
[y3, x3] = ksdensity(ensstates1_1(varplot_1,:,timeplot_1));
[y4, x4] = ksdensity(statesA1_1(varplot_1,:,timeplot_1));

%now for benchmark - mehtod b2
[v1, u1] = ksdensity(tempensstates3_1(varplot_1,:,timeplot_1));
[v2, u2] = ksdensity(eysample3_1(varplot_1,:,timeplot_1-1));  %NOTE THEY DON'T ALIGN PERFECTLY, HENCE NEED TO HAVE -1.
[v3, u3] = ksdensity(ensstates3_1(varplot_1,:,timeplot_1));
[v4, u4] = ksdensity(statesA3_1(varplot_1,:,timeplot_1));

%now for second benchmark - method b1 
[t1, s1] = ksdensity(tempensstates2_1(varplot_1,:,timeplot_1));
[t2, s2] = ksdensity(adderr2_1(varplot_1,:,timeplot_1));  
[t3, s3] = ksdensity(ensstates2_1(varplot_1,:,timeplot_1));
[t4, s4] = ksdensity(statesA2_1(varplot_1,:,timeplot_1));

%determine xlims
xbuff = 0.0005
xbuff2 = 0.00005

%CASE STUDY 1:
xlims1 = round([min([min(x1), min(u1),min(s1), x_true_1(varplot_1,timeplot_1)]) - xbuff max([max(x1), max(u1), max(s1), x_true_1(varplot_1,timeplot_1)]) + xbuff],2);
xlims2 = round([min([min(x2), min(u2), min(s2)]) max([max(x2), max(u2), max(s2), min(s2)])],2);
xlims3 = round([min([min(x3), min(u3), min(s3), x_true_1(varplot_1,timeplot_1)])- xbuff max([max(x3), max(u3), max(s3), x_true_1(varplot_1,timeplot_1)])+ xbuff],2);
xlims4 = round([min([min(x4), min(u4), min(s4), x_true_1(varplot_1,timeplot_1)]) max([max(x4), max(u4), max(s4), x_true_1(varplot_1,timeplot_1)])],3);

%MANUAL LIMS:
%xlims1 = ([2.35 2.5]);
xlims2 = ([-0.08 0.08]);
xlims3 = ([-2.02 -1.9]);

figure
subplot(norows,nocols,1)
h2 = plot(x1, y1, 'r', u1, v1, 'k', s1, t1, 'b', x_true_1(varplot_1, timeplot_1)*ones(1, 2), [0 max([max(y1), max(v1), max(t1)])], '--g');
%now change colour of green line:
set(h2(4), 'Color', grycol)
set(h2, 'LineWidth', 2)
hold on 
ylims = ylim;
h4 = plot(x_true_1(varplot_1,timeplot_1)*ones(1,2), ylims ,'--g', 'LineWidth', 2);
set(h4, 'Color', grycol)
xlim([xlims1(1) - xbuff, xlims1(2) + xbuff])
%xlabel(['x', num2str(varplot_1), '(t)'])
xlabel('{\it{M}}({\bf{x}}_{\it{j-1}}^{\it{a}}[k])')
%ylabel(['p(x', num2str(varplot_1), '(t))'])
ylabel('{\it{p}}({\it{M}}({\bf{x}}_{\it{j-1}}^{\it{a}}[k]))')
title('Modelled')
set(gca, 'XTick', xlims1(1):0.02:xlims1(2))
grid on 
set(gca, 'FontSize', fs, 'Position',[0.13 0.194560669456067 0.213405797101449 0.730439330543931],'TickLength',[0.025 0.025], 'LineWidth', 1.8,'TitleFontSizeMultiplier', 0.9);


subplot(norows,nocols,2)
h3 = plot(x2, y2, 'r', u2, v2, 'k', s2, t2, 'b');
set(h3, 'LineWidth', 2)
xlim([xlims2(1) - xbuff, xlims2(2) + xbuff])
xlabel('{\bf{\eta}}_{\itj}[k]')
ylabel('{\itp}({\bf{\eta}}_{\itj}[k])')
title('Error')
set(gca, 'XTick', xlims2(1):0.04:xlims2(2))
grid on
set(gca, 'FontSize', fs, 'Position',[0.380592139205477 0.194560669456067 0.213405797101449 0.728347280334727],'TickLength',[0.025 0.025], 'LineWidth', 1.8,'TitleFontSizeMultiplier', 0.9);


subplot(norows,nocols,3)
h4 = plot(x3, y3, 'r', u3, v3, 'k', s3, t3, 'b', x_true_1(varplot_1, timeplot_1)*ones(1, 2), [0 max([max(y3), max(v3), max(t3)])], '--g');
%now change colour of green line:
set(h4(4), 'Color', grycol)
set(h4, 'LineWidth', 2)
legend([legnames, 'True'], 'Location', 'south', 'Orientation', 'Horizontal', 'FontSize', 24)
xlim([xlims3(1) - xbuff, xlims3(2) + xbuff])
hold on 
ylims = ylim;
h5 = plot(x_true_1(varplot_1,timeplot_1)*ones(1,2), ylims ,'--g', 'LineWidth', 2)
set(h5, 'Color', grycol)
xlabel('{\bf{x}}_{\itj}^{\itf}[k]')
ylabel('{\it{p}}({\bf{x}}_{\itj}^{\itf}[k])')
title('Forecast')
set(gca, 'XTick', xlims3(1):0.04:xlims3(2))
grid on 
set(gca, 'FontSize', fs, 'Position',[0.629566155433608 0.186192468619247 0.213405797101449 0.734623430962342],'TickLength',[0.025 0.025], 'LineWidth', 1.8,'TitleFontSizeMultiplier', 0.9);


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





