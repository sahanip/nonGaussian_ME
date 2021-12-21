function [diff, ey_est, xcovs, eyall,xbinc,xbinind,xbins, mest] = L96testFSOLVEF_window2021_withbin(xi, exhid, F, K, tstep, yobs, obsind,hidind, windl, obsfreq, bw)


%xi= K x 1  matrix of states at time t - only accepts one state, can't have
%ensemble - best possible IC 
%F, K, tstep = parameters for Lorenz
%yobs = n x windl matrix of observed variables for time window (including initial time)
%windl = no. of time steps in advance to consider 
%exhid = m x windl matrix estimate of errors on hidden states. 
%bw = bandwidth of the gaussian kernel i.e. standard deviation

%STEP 1: first calculate full state and error vector based on ic (xi) and error on
%hidden states (exhid):
xstart= NaN*ones(K,windl+1);
xstart(:,1) = xi;

for k = 1:windl
    x0 = xstart(:,k);
    x0(obsind) = yobs(:,k);
    %simulate to next obs:
    xtemp = NaN*ones(K, obsfreq+1);  xtemp(:,1) = x0;
    for m = 1:obsfreq 
        xtemp(:,m+1) = lorenz96(xtemp(:,m),F,K, tstep,1);
    end
    xstart(:,k+1) = xtemp(:,end);
    xstart(hidind,k+1) = xstart(hidind,k+1) + exhid(:,k);
end

yi1 = xstart(obsind,2:end);

%calculate euclidean distance (l2 norm):
ey_est = yobs(:,2:end) - yi1;

%OPTION3: assuming cyclic symmetry and state dependence:
eyall = NaN*ones(K, windl);
eyall(obsind,:) = ey_est;
eyall(hidind,:) = exhid;
xcovs = xstart(:,1:end-1);

Xd = xcovs(:);
Yd = eyall(:);
diff = 0;
%now estimate conditional variance using kernel:
mest = Xd*0;  %preallocate

for i = 1:length(Xd)
    % calc cond mean:
    K = normpdf(Xd, Xd(i), bw);
    normK = K/sum(K);
    %for j = 1:length(Xd)
        %diff = diff + ((Yd(i) - Yd(j))^2)*normpdf(Xd(j), Xd(i), bw);
        %cmean = cmean + normpdf(Xd(j), Xd(i), bw)*Yd(j);
    %end
    diff = diff + (Yd(i) - sum(normK.*Yd))^2;
    mest(i) = sum(normK.*Yd);
end
diff = (1/length(Xd))*diff;


%following is just for plotting figure
%now bin things - bin such that there are at least 5 in each bin:
xbins = [floor(min(min(xcovs))):ceil(max(max(xcovs)))];  
[xbinc, xbinind] = histc(xcovs(:), xbins);

%xbinc = count in each bin 
%xbinind = separator for each bin.  

%now check if sufficient no. in each bin, if not, aggregate:
diff = 0;
for j = 1:length(xbins)-1
    if xbinc(j) > 5        
        diff = diff + var(eyall(find(xbinind == j)));
    else
        %add to next bin:
        xbinind(find(xbinind == j)) = j+1;
        xbinc(j+1) = xbinc(j+1) + xbinc(j); xbinc(j) = 0;
    end
end        





