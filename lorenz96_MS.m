%multi-scale L96, not periodicity on y is across full range, as per
%Crommelin and vanden-eijnden

function rhs=lorenz96_MS(t,y,F, hx, hy, eps_g, J,K) 

%y = vector of size (J+1)*K for all X and Y vars. Grouping: y(1) = X1,
%y(2:J) = Y1,1:YJ,1 (i.e. all y's associated to J) and so on. 
%K = no slow vars
%J = no fast vars
%eps_g = epsilon 

counter = 1;
for m = 1:K
    xs(m,1) = y(counter);
    ys(:,m) = y(counter+1:counter+J);
    counter = counter+J+1;
end

rhs_x =[xs(end);xs(1:end-1)].*([xs(2:end);xs(1)] - [xs(end-1:end);xs(1:end-2)]) - xs + F*xs.^0 + (hx/J)*(sum(ys)');

ysnew = ys(:);
xsnew = repmat(xs', J,1);
xsnew = xsnew(:);

rhs_y = (1/eps_g)*(([ysnew(2:end); ysnew(1)]).*([ysnew(end); ysnew(1:end-1)] - [ysnew(3:end); ysnew(1:2)]) - ysnew + hy*xsnew);

%now combine into original form:
rhs_y2 = reshape(rhs_y, J,K);
rhs = [rhs_x'; rhs_y2];
rhs = rhs(:);




