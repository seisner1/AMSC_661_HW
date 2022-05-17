



%---------problem 2: Gudonov--------------------


dt = 0.001;
N = 100;
h = 4./N;

rho = zeros(1,N);

T = 12;

i = 1:100;

x = 0 + i.*h;

rho(x <= 1) = 1;

% rho_L = rho(end);
% rho_R = rho(1);
% 
% rho = [rho_L,rho,rho_R];

t = 0;
close all
figure;
hold on;
R = mod(i,N)+1;
L = mod(i-2,N)+1;
while t <= T
    
    rho(i) = rho(i) - (dt./h).*(f_u_star(rho(i),rho(R)) - f_u_star(rho(L),rho(R)));
    
%     rho_L = rho(end);
%     rho_R = rho(1);
% 
%     rho = [rho_L,rho(i),rho_R];
    
    t = t + dt;
    
    u_exact = (x + rho(i)*t)/t;
    
    if abs(t - 7) <= 1e-6
        figure(1);
        grid;
        xlabel('x');
        ylabel('u(x)');
        title('u(x,t=7)');
        plot(x,rho(i));
        plot(x,u_exact);
        
    end
    if abs(t - 11) <= 1e-6
        figure(2);
        hold on;
        xlabel('x');
        ylabel('u(x)');
        title('u(x,t=11)');
        plot(x,rho(i));
        plot(x,u_exact);
    end
    
%     if abs(t - 10) <= 1e-6
%         figure(3);
%         plot(x,rho(i));
%         xlabel('x');
%         ylabel('rho');
%         title('rho(x,10)');
%     end
%     plot(x,rho(i));
%     drawnow;
end

spectral();

function out = f(rho)

    out = (1/2)*rho.^2;


end


function out = f_u_star(rho_L,rho_R)
    
    
    entropy = (f(rho_L) - f(rho_R))./(rho_L - rho_R);
    
    flag = entropy > zeros(1,length(entropy));
    
    rho_star = flag.*rho_L + (~flag).*rho_R;
    
    out = f(rho_star);


end

%problem 3: spectral solver

function spectral()
% solves u_t + u_{xxx} + (0.5u^2)_x = 0, i.e.,
% u_t = -u_{xxx} - (0.5u^2)_x

init_data = 2;

N = 512;
L = 4;
x = linspace(0,L,N+1);
x(N + 1) = [];
k = -N/2 : (N/2 - 1);
u = zeros(1,N);
% initial data
if init_data == 1 
    u0 = u;
end
if init_data == 2
    u0 = u;
    u0(x <= 1) = 1;
end
dt = 0.1; % time step
% figure; clf; 
% hpic = plot(x,u0,'LineWidth',2,'color','r'); % plot of the numerical solution
% hold on;
grid
% drawnow
if init_data == 1
    % plot of the exact solution for u_0(x) given by  (*)
    hp = plot(x,u0,'LineWidth',2); 
    axis([-L/2 L/2 -0.01 1.01]);
end
%
tmax = 12;
t = 0;
freq = k.*(2*pi/L); % frequencies
freq3 = freq.^3;
freq2 = freq.^2;
freq4 = freq.^4;
e3=exp((-0.01.*freq2)*dt);

U = u;
while (t<tmax) 
    t=t+dt;
    vhat=fftshift(fft(u0)); % v in the Fourier space
    k1=rhs(0,vhat);
    k2=rhs(0.5*dt,vhat+0.5*dt*k1);
    k3=rhs(0.5*dt,vhat+0.5*dt*k2);
    k4=rhs(dt,vhat+dt*k3);
    vhat_new=vhat+dt*(k1+2*k2+2*k3+k4)/6;
    unew=ifft(ifftshift(e3.*vhat_new)); % return to u in the x-space
%     set(hpic,'xdata',x,'ydata',real(unew));
    if init_data == 1
        y = -N/2 + mod(x - t/3 + N/2,N);
        %set(hp,'xdata',x,'ydata',1./(cosh((y)/sqrt(12))).^2);
        axis([-N/2 N/2 -0.01 1.01]);
    end
    u0=unew;
%     drawnow
    
    if abs(t - 7) <= 1e-6
        figure(1);
        plot(x,real(unew));
        legend('numerical','analytic','spectral');
    end
    if abs(t - 11) <= 1e-6
        figure(2);
        plot(x,real(unew));
        legend('numerical','analytic','spectral');
    end
    
    U(end+1,:) = unew;
end
%     figure;
%     xlabel('x');
%     ylabel('t');
%     imagesc(x,0:0.1:200,real(U));
%     xlabel('x');
%     ylabel('t');
%     set(gca, 'YDir', 'normal')
%     colorbar
end


function RHSvhat=rhs(dt,vhat)
% v should be a row vector
% RHSvhat = - e^{-tL}(1i*k*hat{(e^{tL}v)^2/2} 
N=size(vhat,2);
L = N;
k=-N/2 : (N/2 - 1);
freq =k.*(2*pi/L);
freq2 = freq.^2;
freq4 = freq.^4;
e3=exp((-0.01.*freq2)*dt);
em3=exp(-(-0.01.*freq2)*dt);
    vhat1=vhat.*e3;          % e^{tL}v in the Fourier space 
    v1=ifft(ifftshift(vhat1));      % exp(tL)v in the x-space
    v2=0.5*v1.^2;          % [exp(tL)v]^2 in the x-space
    RHSvhat=-em3.*(1i*freq).*fftshift(fft(v2)); % exp(-tL)[[(exp(tL)v)]_x] in the Fourier space
end

function u0 = initfunc(x)

    u = zeros(size(x));
    
    u(x<=1) = 1;
    u(x<0) = 0;
    
    u0 = u;


end