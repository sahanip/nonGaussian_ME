%compute data for ACF for case study 1 method B2 

clear all
close all

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%USER INPUTS:

load bwresults_CS1_B2.mat

load evaldata_syn_L96_CS1.mat

clearvars -except bw bwm covindex datatrain Fcons eps_g hx hy J K obsfreq ts x0_true y0_true 

ballsize = [1,0.01]; %for xk, xk-1,exk directions.
pdensthresh = 10;

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

tstep = ts;
novars = K;

nodatapts = 10^5;
simlength = nodatapts;

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
    
xtest(:,2) = tempx; eysample(:,2) = x_truefull(:,2) - tempx;  

fixedres =1;
n = 1;

errbnds = [min(datatrain(:,1))-0.05, max(datatrain(:,1))+0.05]; 

maxbwfac = 1.5;
initbw = bw;  %fixed multivariate bandwidth estimate from np (only for calculating the bandwidth factors)
initbwm = bwm;
[~,~, bwfactors] = CKDE_varbw_silv_estbw(datatrain, covindex, initbw, initbwm, 1/2);
bwfactors(find(bwfactors > maxbwfac)) = maxbwfac;
bwadapt = repmat(bwfactors, 1, length(bw)).*repmat(bw, size(datatrain,1), 1);  %using the multivariate bandwidth estimate 

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
        
        arcov = max([eysample(p,i), errbnds(1)]);
        arcov = min([arcov, errbnds(2)]);  
        covaltemp = [xtest(p,i), arcov];
        temperrest = sampleCKDE(datatrain(:,[fixedres,covindex(fixedres,1:size(covindex,2))]), bwadapt, n, covindex, covaltemp, 0,0);

        if isempty(temperrest)                
            %redefine covariates: locate closest point in covariate space and then sample:
            [~,mc] = min(sum((datatrain(:, covindex) - repmat(covaltemp, size(datatrain,1),1)).^2,2));
            covaltemp = datatrain(mc,covindex);
            temperrest = sampleCKDE(datatrain(:,[fixedres,covindex(fixedres,1:size(covindex,2))]), bwadapt, n, covindex, covaltemp, 0,0);
        end        
        eysample(p,i+1) = temperrest;
        xtest(p,i+1) = sttemp(p,end) + eysample(p,i+1); 
        
        %now check if reasonable:
        [xtest(p,i+1), eysample(p,i+1)] = resample_region_v2(i, p, xtest, eysample, datatrain, bwadapt,covindex, fixedres, ballsize,pdensthresh,sttemp(p,end));
       
        clear temperrest       
    end
       
end

%save results:
save('ACFdata_CS1_B2.mat', 'xtest', 'eysample', 'x_truefull') 

%now do quick acf calc:

%now calculate ACF for each 
maxlag = 1000;

for k = 1:K
    acf1(k,:) = autocorr(xtest(k,:), maxlag);
    acftrue(k,:) = autocorr(x_truefull(k,:), maxlag);    
end

figure
plot([0:maxlag], mean(acftrue), 'k', [0:maxlag], mean(acf1), 'r')
legend('True', 'Estimated')



    

