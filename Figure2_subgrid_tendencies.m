%Figure1: Sub-grid tendencies of the 2 different regimes considered in this
%study 

clear all
close all

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%USER INPUTS:

%load data used for case study 1
load traininputdata_syn_L96_CS1.mat
x1_sg1 = x1_sg; y1_sg1 = y1_sg; 
clearvars -except *1

%load data used for case study 2
load traininputdata_syn_L96_CS2.mat
x1_sg2 = x1_sg; y1_sg2 = y1_sg; 
clearvars -except *1 *2

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

figure
subplot(1,2,1)
plot(x1_sg1(:), y1_sg1(:), '.')
xlabel('{\bf{X}}[k]'); ylabel(['Sub grid tendency ', '{\bf{U}}[k]'])
xlim([-10 15])
title('Case Study 1')
set(gca, 'FontSize', 24,'TickLength',[0.025 0.025], 'LineWidth', 1.8,'TitleFontSizeMultiplier', 0.9);


subplot(1,2,2)
plot(x1_sg2(:), y1_sg2(:), '.')
set(gca, 'YTick', [-10:2.5:5])
%xlabel('X[k]'); ylabel('Sub grid tendency U[k]')
xlim([-10 15])
title('Case Study 2')
set(gca, 'FontSize', 24,'TickLength',[0.025 0.025], 'LineWidth', 1.8,'TitleFontSizeMultiplier', 0.9);

annotation(gcf,'textbox',...
    [0.4306066252588 0.81924882629108 0.0310910973084887 0.0985915492957746],...
    'String',{'a)'},...
    'LineStyle','none',...
    'FontSize',24,...
    'FontName','Helvetica',...
    'FitBoxToText','off');

annotation(gcf,'textbox',...
    [0.871600414078678 0.81924882629108 0.0310910973084887 0.0985915492957746],...
    'String','b)',...
    'LineStyle','none',...
    'FontSize',24,...
    'FontName','Helvetica',...
    'FitBoxToText','off');




