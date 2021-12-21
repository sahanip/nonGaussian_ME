function [CRPSval, intsum] = CRPS(yfor, wghts, yobs)

%This function calculates the CRPS, considers weights if the members are
%not equally weighted. 

%INPUTS:
%yfor = n x t matrix of ensemble forecasts for each time 
%wghts = n x t matrix of associated weights on the forecast members
%yobs = 1 x t vector of obs

t = length(yobs);

Ffor = [0.01:0.01:1]';

for i = 1:t
    
    xfor = NaN*ones(length(Ffor), 1);
    
    for j = 1:length(Ffor)
        xfor(j,1) = weightquantile(yfor(:,i)', wghts(:,i)',Ffor(j));
    end
    
   
    %estimate integral:
    %separate computation into 2 parts for the heaviside:
    b = find(xfor >= yobs(i));
    
    if yobs(i) <= xfor(1)
        %yobs less than range of forecast
        intsum(i) = trapz(xfor, (1-Ffor).^2) + (xfor(1) - yobs(i));
        
    elseif yobs(i) >= xfor(end)
        %yobs larger than range of forecast 
        intsum(i) = trapz(xfor, Ffor.^2) + (yobs(i) - xfor(end));
    else       
        
        if isempty(find(xfor == yobs(i)))
            a1 = find(xfor > yobs(i));  a2 = find(xfor < yobs(i));
            c = interp1(xfor([a2(end) a1(1)]), Ffor([a2(end) a1(1)]), yobs(i));   %to avoid issues with non strictly monotonically increasing variables 
            xfor1 = [xfor(1:b(1)-1); yobs(i)];
            Ffor1 = ([Ffor(1:b(1)-1); c]).^2;  %integrand, left of Heaviside step
            xfor2 = [yobs(i); xfor(b(1):end)];
            Ffor2 = (1 - [c; Ffor(b(1):end)]).^2;  %integrand, right of Heaviside step.             
        else            
            xfor1 = xfor(1:b(1));
            Ffor1 = Ffor(1:b(1)).^2;
            xfor2 = xfor(b(1):end);
            Ffor2 = (1 - Ffor(b(1):end)).^2;
        end
        intsum(i) = trapz(xfor1, Ffor1) + trapz(xfor2, Ffor2);
    end
end

CRPSval = mean(intsum);



