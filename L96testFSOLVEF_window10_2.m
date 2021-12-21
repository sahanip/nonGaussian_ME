function [diff,xstart, J_err, J_obs] = L96testFSOLVEF_window10_2(xi, exguess, F, K, tstep, yobs, obsind, windl, obsfreq, RI, QI)

%INPUTS: 
%xi= K x 1  matrix of states at time t - only accepts one state, can't have
%ensemble - best possible IC 
%F, K, tstep = parameters for Lorenz
%yobs = n x windl matrix of observed variables for time window (including initial time)
%windl = no. of time steps in advance to consider 
%exhid = m x windl matrix estimate of errors on hidden states. 
%QI = K x K inverse of process noise covariance matrix 
%RI = n x n inverse of observation error covariance matirx 
%exguess = K x windl matrix of guesses (for observed variables also)

xstart= NaN*ones(K,windl+1);
xstart(:,1) = xi;

for k = 1:windl
    x0 = xstart(:,k);
    %x0(obsind) = yobs(:,k);
    %simulate to next obs:
    xtemp = NaN*ones(K, obsfreq+1);  xtemp(:,1) = x0;
    for m = 1:obsfreq 
        xtemp(:,m+1) = lorenz96(xtemp(:,m),F,K, tstep,1);
    end
    xstart(:,k+1) = xtemp(:,end);
    xstart(:,k+1) = xstart(:,k+1) + exguess(:,k);
end

%OPTION 3: WEAK CONSTRAINT 4DVAR TYPE SETUP:

for m = 1:windl
    J_err(m) = 0.5*(exguess(:,m)'*QI*exguess(:,m));
    J_obs(m) = 0.5*(((xstart(obsind,m+1) - yobs(:,m+1))')*RI*(xstart(obsind,m+1) - yobs(:,m+1))); 
end

diff = sum(J_err + J_obs);





