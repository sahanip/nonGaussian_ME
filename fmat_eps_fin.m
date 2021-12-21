function [kx, ky, Ffluc] = fmat_eps_fin(X0,Y0, F,K, J, eps_g, hx, hy)


kx = zeros(K,1);   % for the slow variables 
%ky = zeros(J*K, 1);  %for the fast variables 

if length(F) == 1
    F = F*ones(size(kx));
end

%1. First the slow variables:
rhs_x =[X0(end);X0(1:end-1)].*([X0(2:end);X0(1)] - [X0(end-1:end);X0(1:end-2)]) - X0;

%add the F and the fast variables  
Y0new = reshape(Y0, J, K);
Ffluc = (hx/J)*(sum(Y0new)');
if or(max(isnan(Ffluc)) ==1, max(isinf(abs(Ffluc))) == 1)
    error('Note, error in Y0!')
end
kx = rhs_x + F + Ffluc;


%2. Now for the fast variables:

xsnew = repmat(X0', J,1);
xsnew = xsnew(:);

ky = (1/eps_g)*(([Y0(2:end); Y0(1)]).*([Y0(end); Y0(1:end-1)] - [Y0(3:end); Y0(1:2)]) - Y0 + hy*xsnew);




