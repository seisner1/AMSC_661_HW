function KdVrkm()
% solves u_t + u_{xxx} + (0.5u^2)_x = 0, i.e.,
% u_t = -u_{xxx} - (0.5u^2)_x

init_data = 2;

close all
N = 512;
L = 32*pi;
x = linspace(0,L,N+1);
x(N + 1) = [];
k = -N/2 : (N/2 - 1);
u = zeros(1,N);
% initial data
if init_data == 1 
    u0 = 1./(cosh(x/sqrt(12))).^2;
end
if init_data == 2
    u0 = cos(x/16).*(1 + sin(x/16));
end
dt = 0.1; % time step
%figure; clf; 
%hpic = plot(x,u0,'LineWidth',2,'color','r'); % plot of the numerical solution
%hold on;
%grid
%drawnow
if init_data == 1
    % plot of the exact solution for u_0(x) given by  (*)
    hp = plot(x,u0,'LineWidth',2); 
    axis([-L/2 L/2 -0.01 1.01]);
end
%
tmax = 200;
t = 0;
freq = k.*(2*pi/L); % frequencies
freq3 = freq.^3;
freq2 = freq.^2;
freq4 = freq.^4;
e3=exp((freq2 - freq4)*dt);

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
    %set(hpic,'xdata',x,'ydata',real(unew));
    if init_data == 1
        y = -N/2 + mod(x - t/3 + N/2,N);
        %set(hp,'xdata',x,'ydata',1./(cosh((y)/sqrt(12))).^2);
        axis([-N/2 N/2 -0.01 1.01]);
    end
    u0=unew;
    %drawnow
    
    U(end+1,:) = unew;
end
    figure;
    xlabel('x');
    ylabel('t');
    imagesc(x,0:0.1:200,real(U));
    xlabel('x');
    ylabel('t');
    set(gca, 'YDir', 'normal')
    colorbar
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
e3=exp((freq2-freq4)*dt);
em3=exp(-(freq2-freq4)*dt);
    vhat1=vhat.*e3;          % e^{tL}v in the Fourier space 
    v1=ifft(ifftshift(vhat1));      % exp(tL)v in the x-space
    v2=0.5*v1.^2;          % [exp(tL)v]^2 in the x-space
    RHSvhat=-em3.*(1i*freq).*fftshift(fft(v2)); % exp(-tL)[[(exp(tL)v)]_x] in the Fourier space
end
