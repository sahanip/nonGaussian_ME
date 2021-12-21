%compute data for ACF case study 2 for proposed approach 

clear all
close all

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%USER INPUTS:

load bwresults_CS2_PM.mat  
clearvars -except bw bwm covindex datatrain 

load evaldata_syn_L96_CS2.mat 

clearvars -except  Fcons eps_g hx hy J K obsfreq ts x0_true y0_true bw bwm covindex datatrain 

%FOLLOWING FOR THE REJECT RESAMPLE:
ballsize = [1,1,0.05]; %for xk, xk-1,exk directions.
pdensthresh = 10;

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

tstep = ts;
novars = K;

nodatapts = 10^5;
simlength = nodatapts

%use initial condition from evaluation period:
x_truefull = NaN*ones(K, simlength); x_truefull(:,1) = x0_true; 
xtest = x_truefull; eysample = NaN*ones(K, simlength);  

%determine initial error:
%determine true error at first time step:
tempx = x_truefull(:,1);
for j = 1:obsfreq
    tempx = lorenz96(tempx,Fcons,K, tstep,1);
end

ytmp2 = [x0_true'; reshape(y0_true, J,K)];
ytmp2 = ytmp2(:);

options = odeset('Maxstep',0.01,'RelTol',1e-5,'AbsTol',1e-4);

%now advance full model:
[ttmp,ytmp]=ode45(@(t,y) lorenz96_MS(t,y,Fcons, hx, hy, eps_g, J,K), [0 obsfreq*ts], ytmp2,options); 
ytmp2 = ytmp(end,:)';  
x_truefull(:,2) = ytmp2(1:(J+1):end);
    
eysample(:,2) = x_truefull(:,2) - tempx;
xtest(:,2) = tempx; 

fixedres =1;
n = 1;
ptsamp = 100;

bwadapt = repmat(bw, size(datatrain,1), 1);  %no adaptive bandwidth 

for i = 2:simlength-1
    i
    
    %generate full model run:    
    [ttmp,ytmp]=ode45(@(t,y) lorenz96_MS(t,y,Fcons, hx, hy, eps_g, J,K), [0 obsfreq*ts], ytmp2,options); 
    ytmp2 = ytmp(end,:)';  
    x_truefull(:,i+1) = ytmp2(1:(J+1):end);
        
    %generate model forecast:
    sttemp = NaN*ones(K, obsfreq+1); 
    sttemp(:,1) = xtest(:,i);  
    
    for j = 1:obsfreq
        sttemp(:,j+1) = lorenz96(sttemp(:,j),Fcons,K, tstep,1);
    end
  
    %now generate additive errors:
        
    for p = 1:novars
        
        if p == 1
            pp = novars;
        else
            pp = p - 1;
        end
        
        %First for my method:
        eysample(p,i+1) = sampleCKDE(datatrain(:,[fixedres,covindex(fixedres,1:size(covindex,2))]), bwadapt, n, covindex, [xtest(p,i), xtest(pp,i), eysample(p,i)], 0,0);
        xtest(p,i+1) = sttemp(p,end) + eysample(p,i+1); 
        
        [xtest(p,i+1), eysample(p,i+1)] = resample_region(i, p, pp, xtest, eysample, datatrain, bwadapt,covindex, fixedres, ballsize,pdensthresh,sttemp(p,end)); 
       
    end

    
end


%save results:
save('ACFdata_CS2_PM.mat', 'xtest', 'eysample', 'x_truefull') 

%now do quick acf calc:

%now calculate ACF for each 
maxlag = 250;

for k = 1:K
    acf1(k,:) = autocorr(xtest(k,:), maxlag);
    acftrue(k,:) = autocorr(x_truefull(k,:), maxlag);    
end

figure
plot([0:maxlag], mean(acftrue), 'k', [0:maxlag], mean(acf1), 'r')
legend('True', 'Estimated')
title(expnameACF)


