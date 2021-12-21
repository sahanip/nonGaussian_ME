function [CLOC] = localisation_full_GG(D,ro_loc)

% Gaspari and Cohn, 1999, QJRMS, 125, 723-757, eq. 4.10

%USER INPUTS:
%ro_loc = radius 
%D = no. vars 

% Determine coefficients along rows of covariance localisation matrix CLOC
dist = [0:(D-1)];
%
% L96 is periodic
dist = min([dist;(D-dist)]);
%

ind1 = find(dist <= ro_loc);
r2 = (dist(ind1) / ro_loc) .^ 2;
r3 = (dist(ind1) / ro_loc) .^ 3;
coeffs(ind1) = 1 + r2 .* (- r3 / 4 + r2 / 2) + r3 * (5 / 8) - r2 * (5 / 3);
ind2 = find(dist > ro_loc & dist <= 2*ro_loc);
r1 = (dist(ind2) / ro_loc);
r2 = (dist(ind2) / ro_loc) .^ 2;
r3 = (dist(ind2) / ro_loc) .^ 3;
coeffs(ind2) = r2 .* (r3 / 12 - r2 / 2) + r3 * (5 / 8) + r2 * (5 / 3)...
             - r1 * 5 + 4 - (2 / 3) ./ r1;

% build covariance localisation matrix         
CLOC = zeros(D, D);
for i = 1 : D
    vcol(i : D) = coeffs(1 : D + 1 - i);
    vcol(1 : i - 1) = coeffs(i : -1 : 2);
    CLOC(i, 1:D) = vcol;
end
