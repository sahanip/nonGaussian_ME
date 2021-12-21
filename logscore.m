function [ls, els] = logscore(x, xpdf, y, zeroval) 

%this function calculates the log score time series.  It assigns a user
%defined value for zero probabilities. 

%INPUTS:
%x = n x t matrix of values at which the predictive density is given, t =
%no. of time points 
%xpdf = n x t matrix of corresponidng pdf values 
%y = 1 x t time series of observations 

t = size(x,2);

for i = 1:t
    %calculate p(y):
    py(i) = interp1(x(:,i), xpdf(:,i), y(i));
    
    if isnan(py(i))
        ls(i) = zeroval;
    else
        ls(i) = log(py(i));
    end

end

els = mean(ls); %expected score function 