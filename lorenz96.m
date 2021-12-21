
function X = lorenz96(X0,F,J,h,n)
% X (txJ) = lorenz4D(tf,forcing,perturbation)
% ensure final time is divisible by h

%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
%INPUTS
%J = no. of variables 
%F = forcing value (usu 8 - exhibits chaotic behaviour)
%X0 = initial condition for state variables, J x n matrix  
%h = time step, standard is 0.05
%n = no. of ensemble members 
%XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

%%% the Lorenz model is: (cyclical)
% dX[j]/dt=(X[j+1]-X[j-2])*X[j-1]-X[j]+F
%J=40;               %the number of variables
%h=0.05;             %the time step

%initialize
X =X0 +rk4mat(X0,h,F,J,n); %solved via RK4





