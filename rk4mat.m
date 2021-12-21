function deltay = rk4mat(X0,h,F,J, n)

%same as rk4 but operates with matrix X
% X[t+1] = rk4(X[t],step)

k1 = fmat(X0,F,J,n);
k2 = fmat(X0+1/2.*h.*k1,F,J,n);
k3 = fmat(X0+1/2.*h.*k2,F,J,n);
k4 = fmat(X0+h.*k3,F,J,n);
deltay= 1/6*h*(k1+2*k2+2*k3+k4);