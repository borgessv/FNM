function [T,X] = rk4(fun,T,X0)
% 4th order Runge-Kutta time integrations scheme
% Author: Vitor Borges Santos - borgessv93@gmail.com
% Date: 22/10/2022

dt = T(2) - T(1);
X = [X0 zeros(size(X0,1),length(T))];
T = T.';
for n = 1:length(T)-1
    k1 = dt*fun(T(n),X(:,n));
    X1 = X(:,n) + 1/2*k1;
    k2 = dt*fun(T(n),X1);
    X2 = X(:,n) + 1/2*k2;
    k3 = dt*fun(T(n),X2);
    X3 = X(:,n) + k3;
    k4 = dt*fun(T(n),X3);

    X(:,n+1) = X(:,n) + 1/6*(k1 + 2*k2 + 2*k3 + k4);
end

end