
%---------- Problem 1.c----------

%----- github link: 

N = 100;
close all

format long;

tspan = [0 1.2];

h = 2/(N+1);


u = zeros(N+2,1);
% set IC

x_L0 = -1;
x_R0 = 1;

x_f0 = 0.5*(x_R0 - x_L0);

x_0_0 = 0.5*(x_L0 + x_R0);

j = 1:N;

eta = -1 + j*h;

eta = eta';

%here use the 2 diff IC's

u0 = 1 - (eta*x_f0 + x_0_0).^2;

y0 = [u0;x_L0;x_R0];



%ODE solver

opts = odeset('RelTol', 1e-15, 'AbsTol', 1e-15);
sol = ode15s(@(t,y) bouss(t,y,N), tspan, y0);

times = 0.1:0.1:1.2;

full_eta_t = deval(sol,times);

u_eta_t = full_eta_t(1:N,:);

x_L_t = full_eta_t(N+1,:);

x_R_t = full_eta_t(N+2,:);

x_f_t = 0.5*(x_R_t - x_L_t);

x_0_t = 0.5*(x_L_t + x_R_t);

figure;
hold on;
j = 1:N;
eta = -1 + j*h;
eta = eta';
xlabel('\zeta');
ylabel('u(\zeta,t)/u_{max}(t)');
title('Evolution of u(\zeta,t)/u_{max}(t) over t = [0.1,1.2]');
for i = 1:length(times)
    
    plot(eta,u_eta_t(:,i)./max(u_eta_t(:,i)));



%odefun
end

plot(eta,1-eta.^2,'k', 'LineWidth', 3);

figure;
hold on;
xlabel('x');
ylabel('u(x,t)');
title('Evolution of u(x,t) over t = [0.1,1.2]');

for i = 1:length(times)
    
    plot(eta*x_f_t(i) + x_0_t(i),u_eta_t(:,i));



%odefun
end



%2nd IC
x_L0 = -1;
x_R0 = 1;

x_f0 = 0.5*(x_R0 - x_L0);

x_0_0 = 0.5*(x_L0 + x_R0);

j = 1:N;

eta = -1 + j*h;

eta = eta';

u0 = 1 - 0.99*cos(2*pi*(eta*x_f0 + x_0_0));

y0 = [u0;x_L0;x_R0];



%ODE solver

opts = odeset('RelTol', 1e-15, 'AbsTol', 1e-15);
sol = ode15s(@(t,y) bouss(t,y,N), tspan, y0);

times = 0.1:0.1:1.2;

full_eta_t = deval(sol,times);

u_eta_t = full_eta_t(1:N,:);

x_L_t = full_eta_t(N+1,:);

x_R_t = full_eta_t(N+2,:);

x_f_t = 0.5*(x_R_t - x_L_t);

x_0_t = 0.5*(x_L_t + x_R_t);

figure;
hold on;
j = 1:N;
eta = -1 + j*h;
eta = eta';
xlabel('\zeta');
ylabel('u(\zeta,t)/u_{max}(t)');
title('Evolution of u(\zeta,t)/u_{max}(t) over t = [0.1,1.2]');
for i = 1:length(times)
    
    plot(eta,u_eta_t(:,i)./max(u_eta_t(:,i)));



%odefun
end

plot(eta,1-eta.^2,'k', 'LineWidth', 3);

figure;
hold on;
xlabel('x');
ylabel('u(x,t)');
title('Evolution of u(x,t) over t = [0.1,1.2]');

for i = 1:length(times)
    
    plot(eta*x_f_t(i) + x_0_t(i),u_eta_t(:,i));



%odefun
end

function dy = bouss(t,y,N)
    
    h = 2/(N+1);
    j = 2:(N+1);
    x_R = y(end);
    x_L = y(end-1);
    u_l = 0;
    u_r = 0;
    u_ext = [u_l;y(1:N);u_r];
    u = y(1:N);
    x_f = 0.5*(x_R - x_L);
    eta = -1 + (1:N)*h;
    eta = eta';
    dudeta_end = (1/(2*h))*(4*y(1) - y(2));
    dudeta_1 = (1/(2*h))*(-4*y(N) + y(N-1));
    dudeta = (1./(2*h)).*(u_ext(j+1) - u_ext(j-1));
    dudetadeta = (1./(h^2)).*(u_ext(j+1) + u_ext(j-1) - 2.*u_ext(j));
    
    du = (1./(x_f.^2)).*(-0.5.*((1 + eta).*dudeta_1 + (1 - eta).*dudeta_end).*dudeta + u.*dudetadeta + dudeta.^2);
    
    dx_l = -dudeta_end./x_f;
    dx_r = -dudeta_1./x_f;
    
    dy = [du;dx_l;dx_r];


end