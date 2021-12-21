function [deltax, deltay, kx1, ky1, kx2, ky2, kx3, ky3, kx4, ky4, Fflucall] = rk4mat_eps_fin(X0,Y0, ts,F,K, J, eps_g, hx, hy)

%same as rk4 but operates with matrix X
% X[t+1] = rk4(X[t],step)

%define increments at each variable for the slow and fast variables:
[kx1, ky1,Ffluc1] = fmat_eps_fin(X0,Y0, F,K, J, eps_g, hx, hy);

[kx2, ky2,Ffluc2] = fmat_eps_fin(X0+1/2*ts*kx1,Y0+1/2*ts*ky1, F,K, J, eps_g, hx, hy);

[kx3, ky3,Ffluc3] = fmat_eps_fin(X0+1/2*ts*kx2,Y0+1/2*ts*ky2, F,K, J, eps_g, hx, hy);

[kx4, ky4,Ffluc4] = fmat_eps_fin(X0+ts*kx3,Y0+ts*ky3, F,K, J, eps_g, hx, hy);

deltax= 1/6*ts*(kx1+2*kx2+2*kx3+kx4);

deltay = 1/6*ts*(ky1+2*ky2+2*ky3+ky4);

Fflucall = [Ffluc1, Ffluc2, Ffluc3, Ffluc4];