
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%NOTE THIS DOES NOT ALLOW ENSEMBLES.  ONLY SINGLE RUN!
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

function [X, Y, kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4, Fflucall] = lorenz96_eps_fin(X0,Y0, F,K,J, ts,eps_g, hx, hy)
% X (txJ) = lorenz4D(tf,forcing,perturbation)
% ensure final time is divisible by h

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%INPUTS
%K = no. of slow variables
%J = no. of fast variables for each slow variable 
%F = forcing value (usu 8 - exhibits chaotic behaviour).  can be a vector 
%X0 = initial condition for slow state variables, K x 1 vector 
%Y0 = initial condition for fast variables,  J*K x 1 vector
%ts = time step, standard is 0.05
%eps_g,b,hx, hy are constants
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

%%% the Lorenz model is: (cyclical)
% dX[j]/dt=(X[j+1]-X[j-2])*X[j-1]-X[j]+F
%J=40;               %the number of variables
%h=0.05;             %the time step

%initialize

[deltax, deltay, kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4, Fflucall] = rk4mat_eps_fin(X0,Y0, ts,F,K, J, eps_g, hx, hy);  %solved via RK4  

X = X0 + deltax;
Y = Y0 + deltay;

