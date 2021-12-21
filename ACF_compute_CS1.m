
%loads ACF Data for Case Study 1 and creates ACF,CCF

clear all
close all

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%USER INPUTS:
%proposed method:
load ACFdata_CS1_PM.mat 
xtest1 = xtest; eysample1 = eysample; expname1 = expnameACF;
clearvars -except *1

%method b1
load ACFdata_CS1_B1.mat 
xtest2 = xtest; eysample2 = eysample; expname2 = expnameACF;
clearvars -except *1 *2

%method b2
load ACFdata_CS1_B2.mat 
xtest3 = xtest; eysample3 = eysample; expname3 = expnameACF;
clearvars -except *1 *2 *3 x_truefull

maxlag = 1000;

savefilename = 'ACFdata_CS1.mat';

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

%now calculate ACF for each 

K = size(x_truefull,1);

for k = 1:K
    acf1(k,:) = autocorr(xtest1(k,:), maxlag);
    acf2(k,:) = autocorr(xtest2(k,:), maxlag);
    acf3(k,:) = autocorr(xtest3(k,:), maxlag);
    acftrue(k,:) = autocorr(x_truefull(k,:), maxlag);    
end

%calculate CCF also:
ccftrue = NaN*ones(1, maxlag+1); ccf1 = ccftrue; ccf2 = ccftrue; ccf3 = ccftrue;

%at zero lag:
ccftrue(1) = corr(x_truefull(1,:)', x_truefull(2,:)');  %zeroth time
ccf1(1) = corr(xtest1(1,:)', xtest1(2,:)');  %zeroth time
ccf2(1) = corr(xtest2(1,:)', xtest2(2,:)');  %zeroth time
ccf3(1) = corr(xtest3(1,:)', xtest3(2,:)');  %zeroth time

for t = 1:(maxlag)
    t
    ccftrue(t+1) = corr(x_truefull(1,1:end-(t-1))', x_truefull(2,t:end)');    
    ccf1(t+1) = corr(xtest1(1,1:end-(t-1))', xtest1(2,t:end)');    
    ccf2(t+1) = corr(xtest2(1,1:end-(t-1))', xtest2(2,t:end)');    
    ccf3(t+1) = corr(xtest3(1,1:end-(t-1))', xtest3(2,t:end)');
end

%now save data:
save(savefilename, 'acf*', 'ccf*', 'xtest*', 'x_truefull', 'expname*', 'maxlag')




