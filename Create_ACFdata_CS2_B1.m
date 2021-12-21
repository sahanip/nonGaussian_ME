%compute data for ACF, case study 2, method b1

clear all
close all

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%USER INPUTS:

load trainresults_CS2_B1_lambda=2.mat  %B AND PVEC DATA
alpha = 0.8;  %tuning parameter for the added model error! 

load evaldata_syn_L96_CS2.mat

clearvars -except alpha bmvec driftall Pmmat xt1all Fcons eps_g hx hy J K obsfreq ts x0_true y0_true 

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

tstep = ts;
novars = K;

nodatapts = 10^5;
simlength = nodatapts;

%use initial condition from evaluation period:
x_truefull = NaN*ones(K, simlength); x_truefull(:,1) = x0_true; 
xtest = x_truefull; eysample = NaN*ones(K, simlength); 

ytmp2 = [x0_true'; reshape(y0_true, J,K)];
ytmp2 = ytmp2(:);

options = odeset('Maxstep',0.01,'RelTol',1e-5,'AbsTol',1e-4);

fixedres =1;
n = 1;

for i = 1:simlength-1
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
    eysample(:,i+1) = -alpha*mvnrnd(bmvec', Pmmat)'; 
    xtest(:,i+1) = sttemp(:,end) + eysample(:,i+1);
       
end

%save results:
save('ACFdata_CS2_B1.mat', 'xtest', 'eysample', 'x_truefull') 

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




    

