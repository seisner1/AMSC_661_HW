
dt = 0.001;
N = 100;
h = 50./N;

rho = zeros(1,N);

T = 10;

i = 2:101;

x = -25 + i.*h;

rho = 0.5 + (0.9./pi).*atan(x);

rho_L = 0.5 + (0.9./pi).*atan(-25);
rho_R = 0.5 + (0.9./pi).*atan(25 + h./2);

rho = [rho_L,rho,rho_R];

t = 0;
close all
while t <= T
    
    rho(i) = rho(i) - (dt./h).*(f_u_star(rho(i),rho(i+1)) - f_u_star(rho(i-1),rho(i+1)));
    
    
    t = t + dt;
    
    if abs(t - 1) <= 1e-6
        figure(1);
        plot(x,rho(i));
        xlabel('x');
        ylabel('rho');
        title('rho(x,1)');
    end
    if abs(t - 5) <= 1e-6
        figure(2);
        plot(x,rho(i));
        xlabel('x');
        ylabel('rho');
        title('rho(x,5)');
    end
    if abs(t - 10) <= 1e-6
        figure(3);
        plot(x,rho(i));
        xlabel('x');
        ylabel('rho');
        title('rho(x,10)');
    end
    
end


function out = f(rho)

    out = -rho.*log(rho);


end


function out = f_u_star(rho_L,rho_R)
    
    
    entropy = (f(rho_L) - f(rho_R))./(rho_L - rho_R);
    
    flag = entropy > zeros(1,length(entropy));
    
    rho_star = flag.*rho_L + (~flag).*rho_R;
    
    out = f(rho_star);


end