%% Fundamentals of Numerical Methods
%
% Prof.: Dr. E.T.A. van der Weide
% Student: Vítor Borges Santos
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Assignment 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equation to be solved:
% d^4u/dx^4 = sin(2*pi*x), 0<=x=<1
%
% Boundary conditions: 
% u(x=0)=u(x=1)=0, 
% du/dx(x=0)=du/dx(x=1)=0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Exact Solution
clear; close all; clc

syms u(x) x
d4u(x) = diff(u,x,4);
du(x) = diff(u,x,1);
f(x) = sin(2*pi*x);
BC = [u(0) == 0, u(1) == 0, du(0) == 0, du(1) == 0];

sol_exact(x) = (dsolve(d4u == f,BC));

x_vec = 0:0.01:1;
fig = figure;
plot(x_vec,sol_exact(x_vec),'--k','linewidth',1.25)
hold on


%% Finite Differences Method
N_vec = [17 33 65 129]; % Number of grid points to be used in FDM 
x_limits = [0 1]; % Actual domain of the differential equation
n_ghost = 1; % Number of ghost points to solve the BCs (if needed)

% Approximation of the derivatives needed using FDM (symbolic):
d4u_order = 4;
d4u_neighbors = 4;
du_order = 1;
du_neighbors = 2;
[d4u,d4u_error] = central_FDM(d4u_order,d4u_neighbors);
[du,du_error] = central_FDM(du_order,du_neighbors);
fprintf('Derivative of order %.i using FDM with %.i neighbors: %s + O(%s) \n', d4u_order,d4u_neighbors, char(d4u), char(d4u_error))
fprintf('Derivative of order %.i using FDM with %.i neighbors: %s + O(%s) \n', du_order,du_neighbors, char(du), char(du_error))
fprintf(['\nNote: In the expressions above, the indexes of "u" do not have relation with the actual discretization points, ' ...
    'they are just illustrative. \nFor example, in a 2nd order approximation of a 1st order derivative around an arbitrary point j, ' ...
    'u1 = u(j-1) and u3 = u(j+1), \nwhile u2 = u(j) does not appear in such approximation.\n\n'])   

% Solution of the differential equation using the obtained approximations:
lgd_entry = cell(length(N_vec),1);
L2 = zeros(length(N_vec),1);
h = zeros(length(N_vec),1);
k_aux = 0;
custom_color = hsv(length(N_vec));
for k = 1:length(N_vec)
    N = N_vec(k); % Number of grid points
    N_eff = N + 2*n_ghost; % Effective number of points
    h(k) = range(x_limits)/(N-1); % Grid step size
    x_vec_FDM = x_limits(1):h(k):x_limits(2); % Domain discretization for FDM 
    syms u_var [1 N_eff]
    eq = sym(zeros(1,N_eff));

    % Algebraic equations for boundary points:
    % See that u_0 = u(1 + n_ghost) and u_N = u(N_eff - n_ghost)
    eq(1) = central_FDM_subs(u_var,1+n_ghost,du,h(k)) == 0; % du/dx(x=0) = 0
    eq(2) = u_var(1+n_ghost) == 0; % u(x=0) = 0
    eq(N_eff-n_ghost) = u_var(N_eff-n_ghost) == 0; % u(x=1) = 0
    eq(N_eff) = central_FDM_subs(u_var,N_eff-n_ghost,du,h(k)) == 0; % du/dx(x=1) = 0

    % Algebraic equations for interior points:
    for i = 2+n_ghost:N_eff-n_ghost-1
        eq(i) = central_FDM_subs(u_var,i,d4u,h(k)) == f(h(k)*(i-(2+n_ghost)));
    end

    sol = struct2cell(solve(eq,u_var));
    sol = double(vpa(sol));
    sol_FDM = sol(1+n_ghost:N_eff-n_ghost); % Removing ghost points from solution
    if k >= length(N_vec)-2
        k_aux = k_aux + 1;
        uh(:,k_aux) = sol_FDM(1:2^(k_aux-1):end); % Auxiliary variable be used in Richardson extrapolation
    end

    % Plotting results:
    plot(x_vec_FDM,sol_FDM,'linewidth',1.25,'color', custom_color(k,:))
    lgd_entry{k} = ['FDM - ',num2str(N_vec(k)),' points'];

    % L2 norm of error:
    L2(k) = sqrt(1/(N-1)*sum((double(sol_exact(x_vec_FDM).') - sol_FDM).^2));
end

fig2 = figure;
loglog(h,L2,'--ob','MarkerFaceColor','b')
grid on
xlabel('$h$','Interpreter','latex','fontsize',14)
ylabel('$L_2$','Interpreter','latex','fontsize',14)


%% Richardson Extrapolation
uh = uh(2:end-1,:);
r = h(2)/h(1);
p = log((uh(:,3) - uh(:,2))./(uh(:,2) - uh(:,1)))/log(r);
C = (uh(:,3) - uh(:,2))./(h(length(N_vec)-2).^p.*(r.^p).*(r.^p - 1));
u_ex = uh(:,3) + r.^p./(1-r.^p).*(uh(:,3) - uh(:,2));

figure(fig)
plot(0:h(length(N_vec)-2):1,[0;u_ex;0],'+r')
grid on
xlabel('$x$','Interpreter','latex','fontsize',14)
ylabel('$u(x)$','Interpreter','latex','fontsize',14)
lgd = legend('Exact',lgd_entry{:},'"Exact" - Richardson','interpreter','latex','fontsize',11);