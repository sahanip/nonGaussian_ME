%Figure: Anomaly densities for proposed method and benchmark for an
%observed and unobserved variable.

clear all
close all

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%USER INPUTS:
%proposed method: 
load assimilateresults_CS1_PM.mat
ensstatesm = ensstates;
clearvars -except ensstatesm

%method B2 
load assimilateresults_CS1_B2.mat
%xfitb = xfit; normfplotb = normfplot; 
ensstatesb = ensstates;
clearvars -except ensstatesm ensstatesb

%method B1 
load assimilateresults_CS1_B1_alpha=0.8.mat

load figure7colormap.mat 

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX


%% plot normalised histogram 

nobins = 30;
obsvarplot = 9
hidvarplot = 1
plottime = simlength;

ylims11 = [-0.04 0.04];   %obs, proposed 
ylims12 = [-0.04 0.04];   %obs, method B2 
ylims13 = [-0.18 0.14];     %obs, method B1 

ylims21 = [-0.7 0.9];   %obs, proposed 
ylims22 = [-0.7 0.9];   %obs, method B2 
ylims23 = [-1 1.4];     %obs, method B1 

nolabs = 9;

fs=20;

nrows = 3
ncols = 2

histnm = NaN*ones(K, nobins, size(ensstates,3));  histcm = histnm; histnb = histnm;  histn = histnm;  histcb = histcm; histc = histcm;

X = repmat([0.5, xpts(:)', plottime+0.5], size(histcm,2),1); 
for k = 1:K
    for t = 2:size(ensstates,3)
        
        %proposed method 
        clear tmp histed 
        tmp = squeeze(ensstatesm(k,:,t)) - x_true(k, t);
        [histnm(k,:,t),histed] = histcounts(tmp,nobins);
        histcm(k,:,t) = (histed(1:end-1)+histed(2:end))/2;  %gives the centre of the bins    
        maxcount(k,t) = max(histnm(k,:,t))/n;
        
        %now for next one:
        clear tmp histed 
        tmp = squeeze(ensstatesb(k,:,t)) - x_true(k, t);
        [histnb(k,:,t),histed] = histcounts(tmp,nobins);
        histcb(k,:,t) = (histed(1:end-1)+histed(2:end))/2;  %gives the centre of the bins    
        maxcount(k,t) = max(histnb(k,:,t))/n;
        
        %now for next one:
        clear tmp histed 
        tmp = squeeze(ensstates(k,:,t)) - x_true(k, t);
        [histn(k,:,t),histed] = histcounts(tmp,nobins);
        histc(k,:,t) = (histed(1:end-1)+histed(2:end))/2;  %gives the centre of the bins    
        maxcount(k,t) = max(histn(k,:,t))/n;
    end
end
histnm = histnm/n;
histnb = histnb/n;
histn = histn/n;


figure 
subplot(nrows,ncols,1)  %proposed, obs var
Y = reshape(repmat(squeeze(histcm(obsvarplot,:,1:plottime)), 2, 1), size(histcm,2), plottime*2);   
C = reshape(repmat(squeeze(histnm(obsvarplot,:,1:plottime)), 2, 1), size(histcm,2), plottime*2);
h = pcolor(X,Y,C);
set(h, 'EdgeColor', 'none');
caxis(clims)
colormap(CustomColormap1)
ylabel('Proposed Method', 'FontWeight', 'Bold')
ylim(ylims11)
set(gca, 'YTick', linspace(ylims11(1), ylims11(2), nolabs));
hold on 
plot([1:plottime], zeros(1, plottime), 'k')
title('Observed Variable')
set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8)
grid on 
%set axis position:
set(gca, 'Position',[0.136369423585922 0.735691767193136 0.327253397893116 0.215735294117647]);

subplot(nrows,ncols,3)  %benchmark - method b2 , obs var 
Y = reshape(repmat(squeeze(histcb(obsvarplot,:,1:plottime)), 2, 1), size(histcb,2), plottime*2);   
C = reshape(repmat(squeeze(histnb(obsvarplot,:,1:plottime)), 2, 1), size(histcb,2), plottime*2);
h = pcolor(X,Y,C);
set(h, 'EdgeColor', 'none');
caxis(clims)
colormap(CustomColormap1)
ylabel('{\it{\bfx}}^{\itf}[k] - {\it{\bfx}}[k]')
ylim(ylims12)
set(gca, 'YTick', linspace(ylims12(1), ylims12(2), nolabs));
hold on 
plot([1:plottime], zeros(1, plottime), 'k')
grid on 
%title('ETKF with LW WC 4D-VAR ')
set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8)
set(gca, 'Position',[0.136369423585922 0.436059414251959 0.327253397893116 0.215735294117647]);



subplot(nrows,ncols,5)  %benchmark - method b1 , obs var 
Y = reshape(repmat(squeeze(histc(obsvarplot,:,1:plottime)), 2, 1), size(histc,2), plottime*2);   
C = reshape(repmat(squeeze(histn(obsvarplot,:,1:plottime)), 2, 1), size(histc,2), plottime*2);
h = pcolor(X,Y,C);
set(h, 'EdgeColor', 'none');
caxis(clims)
colormap(CustomColormap1)
%colorbar        
xlabel('time step'); %ylabel(['x', num2str(obsvarplot), ' predicted - truth']);
ylabel('B1 Method', 'FontWeight', 'Bold')
ylim(ylims13)
set(gca, 'YTick', linspace(ylims13(1), ylims13(2), nolabs));
hold on 
plot([1:plottime], zeros(1, plottime), 'k')
grid on 
set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8)
set(gca, 'Position',[0.136369423585922 0.136427061310782 0.327253397893116 0.215735294117647]);



subplot(nrows,ncols,2)  %proposed hid var
Y = reshape(repmat(squeeze(histcm(hidvarplot,:,1:plottime)), 2, 1), size(histcm,2), plottime*2);   
C = reshape(repmat(squeeze(histnm(hidvarplot,:,1:plottime)), 2, 1), size(histcm,2), plottime*2);
h = pcolor(X,Y,C);
set(h, 'EdgeColor', 'none');
caxis(clims)
colormap(CustomColormap1)
ylim(ylims21)
set(gca, 'YTick', linspace(ylims21(1), ylims21(2), nolabs));
hold on 
plot([1:plottime], zeros(1, plottime), 'k')
grid on 
title('Hidden Variable')
set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8)
set(gca, 'Position',[0.57671033267683 0.735691767193136 0.327253397893117 0.215735294117647]);

subplot(nrows,ncols,4)  %benchmar, hid var
Y = reshape(repmat(squeeze(histcb(hidvarplot,:,1:plottime)), 2, 1), size(histcb,2), plottime*2);   
C = reshape(repmat(squeeze(histnb(hidvarplot,:,1:plottime)), 2, 1), size(histcb,2), plottime*2);
h = pcolor(X,Y,C);
set(h, 'EdgeColor', 'none');
caxis(clims)
colormap(CustomColormap1)
ylabel('B2 Method', 'FontWeight', 'Bold')
ylim(ylims22)
set(gca, 'YTick', linspace(ylims22(1), ylims22(2), nolabs));
hold on 
plot([1:plottime], zeros(1, plottime), 'k')
grid on 
set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8)
set(gca, 'Position',[0.57671033267683 0.436059414251959 0.327253397893117 0.215735294117647]);



subplot(nrows,ncols,6)  %benchmar, hid var
Y = reshape(repmat(squeeze(histc(hidvarplot,:,1:plottime)), 2, 1), size(histc,2), plottime*2);   
C = reshape(repmat(squeeze(histn(hidvarplot,:,1:plottime)), 2, 1), size(histc,2), plottime*2);
h = pcolor(X,Y,C);
set(h, 'EdgeColor', 'none');
caxis(clims)
colormap(CustomColormap1)
colorbar('Location', 'southoutside')  
ylim(ylims23)
set(gca, 'YTick', linspace(ylims23(1), ylims23(2), nolabs));
hold on 
plot([1:plottime], zeros(1, plottime), 'k')
grid on 
set(gca, 'FontSize', fs,'TickLength',[0.025 0.025], 'LineWidth', 1.8)
set(gca, 'Position',[0.57671033267683 0.135306553911205 0.327253397893117 0.216855801517225]);

%annotations 

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








