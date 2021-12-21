%Figure 4: ACFs for both tests.

clear all
close all

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%USER INPUTS:

%CASE STUDY 1:
load ACFdata_CS1.mat
acf1_1 = acf1; acf2_1 = acf2; acf3_1 = acf3; acftrue_1 = acftrue; ccf1_1 = ccf1; ccf2_1 = ccf2; ccf3_1 = ccf3; ccftrue_1 = ccftrue; maxlag_1 = maxlag; x_truefull_1 = x_truefull; xtest1_1 = xtest1; xtest2_1 = xtest2; xtest3_1 = xtest3;

clearvars -except *1

%CASE STUDY 2:
load ACFdata_CS2.mat
acf1_2 = acf1; acf2_2 = acf2; acf3_2 = acf3; acftrue_2 = acftrue; ccf1_2 = ccf1; ccf2_2 = ccf2; ccf3_2 = ccf3; ccftrue_2 = ccftrue; maxlag_2 = maxlag; x_truefull_2 = x_truefull; xtest1_2 = xtest1; xtest2_2 = xtest2; xtest3_2 = xtest3;

obsspace_1 = 0.02; %in MTU
obsspace_2 = 0.04; %in MTU

vardens_1 = 9;  %variable to consider for calculating marginal density.
vardens_2 = 2;

plotccf = 1;  %1 if plot CCF, 0 if not.

greycol = 0.5*ones(1,3);
lwid = 2;

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

%calculate marginal densities:
[ftrue_1, xtrueplot_1] = ksdensity(x_truefull_1(vardens_1,:));
[ftest1_1, xtest1plot_1] = ksdensity(xtest1_1(vardens_1,:));
[ftest2_1, xtest2plot_1] = ksdensity(xtest2_1(vardens_1,:));
[ftest3_1, xtest3plot_1]= ksdensity(xtest3_1(vardens_1,:));

[ftrue_2, xtrueplot_2] = ksdensity(x_truefull_2(vardens_2,:));
[ftest1_2, xtest1plot_2] = ksdensity(xtest1_2(vardens_2,:));
[ftest2_2, xtest2plot_2] = ksdensity(xtest2_2(vardens_2,:));
[ftest3_2, xtest3plot_2]= ksdensity(xtest3_2(vardens_2,:));

fs = 20;

if plotccf == 1
    nocols = 2;
    norows = 3;
    
    figure
    %CASE STUDY 1
    subplot(norows,nocols,1)
    h1 = plot([0:maxlag_1]*obsspace_1, mean(acftrue_1), '--m', [0:maxlag_1]*obsspace_1, mean(acf1_1), 'r',[0:maxlag_1]*obsspace_1, mean(acf2_1), 'b',[0:maxlag_1]*obsspace_1, mean(acf3_1), 'k');
    set(h1(1), 'Color', greycol)
    set(h1, 'LineWidth', 2)
    xlim([0 10])
    xlabel('time lag') 
    grid on 
    ylabel('ACF')
    set(gca, 'FontSize', fs)
    set(gca, 'YTick', [-0.5:0.25:1],'TickLength',[0.025 0.025], 'LineWidth', 1.8)
    title('Case Study 1')
    
    subplot(norows,nocols,3)
    h1 = plot([0:maxlag_1]*obsspace_1, ccftrue_1, '--m', [0:maxlag_1]*obsspace_1, ccf1_1, 'r',[0:maxlag_1]*obsspace_1, ccf2_1, 'b',[0:maxlag_1]*obsspace_1, ccf3_1, 'k');
    set(h1(1), 'Color', greycol)
    set(h1, 'LineWidth', 2)
    xlim([0 10])
    xlabel('time lag') 
    grid on 
    ylabel('CCF')
    set(gca, 'FontSize', fs)
    ylim([-0.5 0.5])
    set(gca, 'YTick', [-0.5:0.1:0.5],'TickLength',[0.025 0.025], 'LineWidth', 1.8)

    %marginal densities:
    subplot(norows,nocols,5)
    h2 = plot(xtrueplot_1, ftrue_1, '--m', xtest1plot_1, ftest1_1, 'r', xtest2plot_1, ftest2_1, 'b', xtest3plot_1, ftest3_1, 'k');
    set(h2(1), 'Color', greycol)
    set(h2, 'LineWidth', 2)
    xlabel('X_k')
    ylabel('p(X_k)')
    grid on 
    set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8)
    xlim([-15 20])
    
    %CASE STUDY 2
    subplot(norows,nocols,2)
    h1 = plot([0:maxlag_2]*obsspace_2, mean(acftrue_2), '--m', [0:maxlag_2]*obsspace_2, mean(acf1_2), 'r',[0:maxlag_2]*obsspace_2, mean(acf2_2), 'b',[0:maxlag_2]*obsspace_2, mean(acf3_2), 'k');
    set(h1(1), 'Color', greycol)
    set(h1, 'LineWidth', 2)
    xlim([0 10])
    xlabel('time lag') 
    legend(h1, {'True', 'Proposed', 'B1', 'B2'}, 'Orientation', 'Horizontal', 'Location', 'southoutside', 'FontSize', fs) 
    grid on 
    ylabel('ACF')
    set(gca, 'FontSize', fs)
    set(gca, 'YTick', [-0.5:0.25:1],'TickLength',[0.025 0.025], 'LineWidth', 1.8)
    title('Case Study 2')
    
    subplot(norows,nocols,4)
    h1 = plot([0:maxlag_2]*obsspace_2, ccftrue_2, '--m', [0:maxlag_2]*obsspace_2, ccf1_2, 'r',[0:maxlag_2]*obsspace_2, ccf2_2, 'b',[0:maxlag_2]*obsspace_2, ccf3_2, 'k');
    set(h1(1), 'Color', greycol)
    set(h1, 'LineWidth', 2)
    xlim([0 10])
    xlabel('time lag') 
    grid on 
    ylabel('CCF')
    set(gca, 'FontSize', fs)
    ylim([-0.6 0.6])
    set(gca, 'YTick', [-0.6:0.2:0.6],'TickLength',[0.025 0.025], 'LineWidth', 1.8)

    %marginal densities:
    subplot(norows,nocols,6)
    h2 = plot(xtrueplot_2, ftrue_2, '--m', xtest1plot_2, ftest1_2, 'r', xtest2plot_2, ftest2_2, 'b', xtest3plot_2, ftest3_2, 'k');
    set(h2(1), 'Color', greycol)
    set(h2, 'LineWidth', 2)
    xlabel('X_k')
    ylabel('p(X_k)')
    grid on 
    set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8)
    xlim([-15 20])
    
else
    nocols = 2;
    norows = 2;
    
    figure
    subplot(norows,nocols,1)
    h1 = plot([0:maxlag_1]*obsspace_1, mean(acftrue_1), '--m', [0:maxlag_1]*obsspace_1, mean(acf1_1), 'r',[0:maxlag_1]*obsspace_1, mean(acf2_1), 'b',[0:maxlag_1]*obsspace_1, mean(acf3_1), 'k');
    set(h1, 'LineWidth', 2)
    xlim([0 10])
    xlabel('time lag') 
    legend(h1, {'True', 'Proposed', 'B1', 'B2'}, 'Orientation', 'Horizontal', 'Location', 'southoutside', 'FontSize', 14)  
    grid on 
    ylabel('ACF')
    set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8)

    %marginal densities:
    subplot(norows,nocols,2)
    h2 = plot(xtrue_1, ftrue_1, '--m', xtest1plot_1, ftest1_1, 'r', xtest2plot_1, ftest2_1, 'b', xtest3plot_1, ftest3_1, 'k');
    set(h2, 'LineWidth', 2)
    xlabel('X_k')
    ylabel('p(X_k)')
    grid on 
    set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8)
    xlim([-15 20])
    
    %CASE STUDY 2:
    subplot(norows,nocols,3)
    h1 = plot([0:maxlag_2]*obsspace_2, mean(acftrue_2), '--m', [0:maxlag_2]*obsspace_2, mean(acf1_2), 'r',[0:maxlag_2]*obsspace_2, mean(acf2_2), 'b',[0:maxlag_2]*obsspace_2, mean(acf3_2), 'k');
    set(h1, 'LineWidth', 2)
    xlim([0 10])
    xlabel('time lag') 
    legend(h1, {'True', 'Proposed', 'B1', 'B2'}, 'Orientation', 'Horizontal', 'Location', 'southoutside', 'FontSize', 14)  
    grid on 
    ylabel('ACF')
    set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8)

    %marginal densities:
    subplot(norows,nocols,4)
    h2 = plot(xtrue_2, ftrue_2, '--m', xtest1plot_2, ftest1_2, 'r', xtest2plot_2, ftest2_2, 'b', xtest3plot_2, ftest3_2, 'k');
    set(h2, 'LineWidth', 2)
    xlabel('X_k')
    ylabel('p(X_k)')
    grid on 
    set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8)
    xlim([-15 20])
end

%figure labels:
annotation(gcf,'textbox',...
    [0.321185614849188 0.846153846153846 0.0210417633410673 0.0725274725274725],...
    'String',{'a)'},...
    'LineStyle','none',...
    'FontSize',fs,...
    'FontName','Helvetica',...
    'FitBoxToText','off');

annotation(gcf,'textbox',...
    [0.601348027842226 0.843956043956044 0.0210417633410674 0.0725274725274725],...
    'String','b)',...
    'LineStyle','none',...
    'FontSize',fs,...
    'FontName','Helvetica',...
    'FitBoxToText','off');

annotation(gcf,'textbox',...
    [0.882090487238975 0.846153846153846 0.0210417633410674 0.0725274725274725],...
    'String','c)',...
    'LineStyle','none',...
    'FontSize',fs,...
    'FontName','Helvetica',...
    'FitBoxToText','off');

annotation(gcf,'textbox',...
    [0.32045086099173 0.402054467271857 0.0210417633410673 0.0725274725274722],...
    'String','d)',...
    'LineStyle','none',...
    'FontSize',fs,...
    'FontName','Helvetica',...
    'FitBoxToText','off');

annotation(gcf,'textbox',...
    [0.601126834540589 0.402054467271857 0.0210417633410672 0.0725274725274723],...
    'String','e)',...
    'LineStyle','none',...
    'FontSize',fs,...
    'FontName','Helvetica',...
    'FitBoxToText','off');

annotation(gcf,'textbox',...
    [0.883272315804363 0.402054467271856 0.0210417633410671 0.0725274725274723],...
    'String','f)',...
    'LineStyle','none',...
    'FontSize',fs,...
    'FontName','Helvetica',...
    'FitBoxToText','off');


annotation(gcf,'textbox',...
    [0.0252486610558531 0.735287420301589 0.0558355649069114 0.0518672199170123],...
    'String',{'  Case','Study 1'},...
    'LineStyle','none',...
    'FontSize',fs,...
    'FontName','Helvetica',...
    'FitBoxToText','off');

annotation(gcf,'textbox',...
    [0.026778882938026 0.290165469082077 0.0558355649069113 0.0518672199170122],...
    'String',{'  Case','Study 2'},...
    'LineStyle','none',...
    'FontSize',fs,...
    'FontName','Helvetica',...
    'FitBoxToText','off');




