% Code to solve the problems of assignment 2 of the Fundamentals of
% Numerical Methods course taken at the University of Twente
%
% Author: Vítor Borges Santos - borgessv93@gmail.com
% Date: 22/10/2022
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Assignment 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Equation to be solved using Finite Element Method and Runge-Kutta:
% d^2u/dx^2 = d^2u/dx^2, 0<=x=<1, t>=0
%
% Boundary conditions: 
% u(0,t) = u(1,0) = 0
% Initial conditions:
% u(x,0) = x^10 - x,
% du/dt(x,0) = 0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc

N_vec = [25 50 100];
L2 = zeros(1,length(N_vec));
for k = 1:length(N_vec)
    N = N_vec(k); % Number of grid spaces
    h = 1/N; % Spatial step-size
    x_vec = 0:h:1; 
    dt = h/8; % Time-step [s]
    T = 0:dt:4;
    T_interest = 0:0.5:4;

    %Initial coditions
    U0 = zeros(N+1,1);
    V0 = zeros(N+1,1);

    % Matrices M and R considering BCs: U(0,t) = U(1,t) = 0
    M = zeros(N+1,N+1);
    M(1,1) = 1;
    M(end,end) = 1;
    R = zeros(N+1,N+1);
    R(1,1) = 1;
    R(end,end) = 1;
    for j = 2:N
        if j >= 2
            M(j,j-1) = h/6;
            R(j,j-1) = 1/h;
        end
        M(j,j) = 4*h/6;
        R(j,j) = -2/h;
        if j <= N
            M(j,j+1) = h/6;
            R(j,j+1) = 1/h;
        end

        U0(j) = ((j-1)*h)^10 - h*(j-1); % Initial U for interior points
    end

    % Time integration of the set of equations:
    X0 = [U0;V0];
    [T,X] = rk4(@(t,X) dynamics(t,X,N,M,R), T, X0);
    U = X(1:N+1,:);
    V = X(N+2:2*(N+1),:);
    
    % L2 norm of error based on the fact U_ex(t=0s) = U_ex(t=4s):
    L2(k) = sqrt(1/(N-1)*sum((U(2:N,T==0) - U(2:N,T==4)).^2));

    % Plotting results:
    figure
    ax = gca;
    hold(ax,'on');
    custom_color = hsv(length(T_interest));
    lgd_entry = cell(length(T_interest),1);
    for i = 1:length(T_interest)
        if i == 1
            plot(x_vec,U(:,T==T_interest(i)),'--k','Linewidth',1.25)
        else
            plot(x_vec,U(:,T==T_interest(i)),'color',custom_color(i,:),'Linewidth',1.25)
        end
        lgd_entry{i} = ['$t=$ ',num2str(T_interest(i),'%.1f'), ' s'];
    end
    grid(ax,'on')
    xlabel(ax,'$x$','Interpreter','latex','fontsize',16)
    ylabel(ax,'$u(x)$','Interpreter','latex','fontsize',16)
    legend(ax,lgd_entry{:},'interpreter','latex','fontsize',13)
    ax.FontSize = 12;
end

figure
loglog(N_vec,L2,'ob','MarkerFaceColor','b')
hold on
N_extrap = [1e0 1e1 1e2 1e3];
L2_extrap = 10.^interp1(log10(N_vec),log10(L2),log10(N_extrap),'linear','extrap');
loglog(N_extrap,L2_extrap,'--r')
grid on
xlabel('$N$','Interpreter','latex','fontsize',16)
ylabel('$L_2$','Interpreter','latex','fontsize',16)
legend('Data','Trend line','interpreter','latex','fontsize',13)
ax = gca;
ax.FontSize = 12;


function Xdot = dynamics(t,X,N,M,R)
Udot = X(N+2:2*(N+1),1);
Vdot = M\(R*X(1:N+1,1));
Xdot = [Udot;Vdot];
end