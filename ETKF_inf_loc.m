function [x_up, K, singmat, x_sim] = ETKF_inf_loc(x_obs, x_sim, R, n,H, tol, lambda1, CLOC)

%based on ETKF function except does mean repulsion additive inflation
%(equiv to multiplicative) and localization 

%THIS IS FOR STATE UPDATING ONLY AND REQUIRES LINEAR OBSERVATION
%OPERATOR!!! 

% %USER INPUTS:
% R = covariance matrix of errors 
% n = ensemble size 
% x_obs = m x 1 vector of obs 
% x_obs_sim = m x n matrix of simulated observed variables 
% x_sim = v x n matrix of simulated variables (all)
% H = m x v obs operator (linear matrix)
% lambda1 = infation factor, if no inflation, use 1.  Note WE TAKE THE
% SQUARE ROOT OF IT! 
% CLOC = mask matrix for localization, if no localization, use ones(v,v)

%-------------------------------------------------------------

%first inflate:
x_sim = repmat(mean(x_sim,2), 1, n) + sqrt(lambda1)*(x_sim - repmat(mean(x_sim,2),1,n));

%STEP 1. ENKF Update without perturbed observations:

D = x_obs - mean(H*x_sim,2);  %observation innovation for mean 

P = cov(x_sim');  %background covariance matrix
%localize matrix:
P = CLOC.*P;

%check invertibility:
V = H*P*H' + R;
[~, S, ~] = svd(V);
sings = max(S);

if or(min(sings)/max(sings) < tol,or(isinf(min(sings)/max(sings)), isnan(min(sings)/max(sings))))
    %can't update         
    singmat = 1;
    x_up = x_sim;
else
    %update
    singmat = 0;
    invmat = inv((V));  %so only have to do the matrix inversion once 
    K = P*H'*invmat; %Kalman gain matrix 
    x_upmean = mean(x_sim,2) + K*D;        
    
%STEP 2. Computate correction to covariance matrix:
    Xf_dash = x_sim - repmat(mean(x_sim,2), 1, n);
    W = (1/(n-1))*Xf_dash'*(H'*inv(R)*H)*Xf_dash;
    %going to use singular value decomposition without removing zero
    %columns 
    [U1,S1,V1] = svd(W);
    T = U1*inv(sqrtm(eye(n) + S1))*U1';
    
    x_up = repmat(x_upmean, 1, n) + Xf_dash*T;
end 






 