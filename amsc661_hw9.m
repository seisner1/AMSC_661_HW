

%--------part 5-----------------

%----github link: 

close all

h = 0.05;
a = sqrt(2);
x = -6:0.05:6;
dt = 0.01*(1/a);

lambda_1 = a;
lambda_2 = -a;

C = [a,-a; 1,1];
Lambda = [a,0;0,-a];
C_inv = 0.5*[1/a,1;-1/a,1];

%----IC
phi = phi_func(x);
psi = zero;

%init solver

u = zeros(1,length(x));
w = zeros(2,length(x));

w(2,:) = gradient(phi,h);

y = C_inv*w;

zeta = y(1,:);
eta = y(2,:);

%---lax-friedrich

i = 1:length(x);

t = 0;

T = 4/a + dt;

x0 = -6;

while t < T
    
    R = mod(i,length(x))+1;
    L = mod(i-2,length(x))+1;
    zeta(i) =  0.5*(zeta(R) + zeta(L)) - lambda_1*(dt/(2*h))*(zeta(R) - zeta(L));
    eta(i) = 0.5*(eta(R) + eta(L)) - lambda_2*(dt/(2*h))*(eta(R) - eta(L));
    
    y = [zeta;eta];
    w = C*y;
    
    j = 2:length(x);
    
    u(1) = 0;
    
    sum = 0;
       
    for j = 2:length(x)
        
        sum = sum + 0.5*h*(w(2,j-1) + w(2,j));
        
        u(j) = sum;
        
    end
    
    
    format long
    
    if (abs(t- 1/(2*a))<=1e-12 || abs(t-1/a)<=1e-12 || abs(t-2/a)<=1e-12 || abs(t-4/a)<=1e-12)
        
        figure;
        hold on
        str = sprintf('Lax-Friedrich, t = %.1f(1/a)',round(t*a,1));
        title(str);
        plot(x,u);
        u_exact = 0.5*(phi_func(x + a*t) + phi_func(x - a*t));
        plot(x,u_exact);
        
        
        
    end
    
    t = t + dt;
    
end





%init solver

u = zeros(1,length(x));
w = zeros(2,length(x));

w(2,:) = gradient(phi,h);

y = C_inv*w;

zeta = y(1,:);
eta = y(2,:);

%---Upwind/Downwind

i = 1:length(x);

t = 0;

T = 4/a + dt;

x0 = -6;

while t < T
    
    R = mod(i,length(x))+1;
    L = mod(i-2,length(x))+1;
    zeta(i) =  -lambda_1*(dt/h)*(zeta(i) - zeta(L)) + zeta(i);
    eta(i) = -lambda_2*(dt/h)*(eta(R) - eta(i)) + eta(i);
    
    y = [zeta;eta];
    w = C*y;
    
    j = 2:length(x);
    
    u(1) = 0;
    
    sum = 0;
       
    for j = 2:length(x)
        
        sum = sum + 0.5*h*(w(2,j-1) + w(2,j));
        
        u(j) = sum;
        
    end
    
    
    format long
    
    if (abs(t- 1/(2*a))<=1e-12 || abs(t-1/a)<=1e-12 || abs(t-2/a)<=1e-12 || abs(t-4/a)<=1e-12)
        
        figure;
        hold on
        str = sprintf('Upwind, t = %.1f(1/a)',round(t*a,1));
        title(str);
        plot(x,u);
        u_exact = 0.5*(phi_func(x + a*t) + phi_func(x - a*t));
        plot(x,u_exact);
        
        
        
    end
    
    t = t + dt;
    
end





%init solver

u = zeros(1,length(x));
w = zeros(2,length(x));

w(2,:) = gradient(phi,h);

y = C_inv*w;

zeta = y(1,:);
eta = y(2,:);

%---Lax-Wendroff

i = 1:length(x);

t = 0;

T = 4/a + dt;

x0 = -6;

while t < T
    
    R = mod(i,length(x))+1;
    L = mod(i-2,length(x))+1;
    zeta(i) =  zeta(i) - lambda_1*(dt/(2*h))*(zeta(R) - zeta(L)) + (lambda_1^2)*(dt^2/(2*h^2))*(zeta(R) - 2*zeta(i) + zeta(L));
    eta(i) = eta(i) - lambda_2*(dt/(2*h))*(eta(R) - eta(L)) + (lambda_2^2)*(dt^2/(2*h^2))*(eta(R) - 2*eta(i) + eta(L));
    
    y = [zeta;eta];
    w = C*y;
    
    j = 2:length(x);
    
    u(1) = 0;
    
    sum = 0;
       
    for j = 2:length(x)
        
        sum = sum + 0.5*h*(w(2,j-1) + w(2,j));
        
        u(j) = sum;
        
    end
    
    
    format long
    
    if (abs(t- 1/(2*a))<=1e-12 || abs(t-1/a)<=1e-12 || abs(t-2/a)<=1e-12 || abs(t-4/a)<=1e-12)
        
        figure;
        hold on
        str = sprintf('Lax-Wendroff, t = %.1f(1/a)',round(t*a,1));
        title(str);
        plot(x,u);
        u_exact = 0.5*(phi_func(x + a*t) + phi_func(x - a*t));
        plot(x,u_exact);
        
        
        
    end
    
    t = t + dt;
    
end






%init solver

u = zeros(1,length(x));
w = zeros(2,length(x));

w(2,:) = gradient(phi,h);

y = C_inv*w;

zeta = y(1,:);
eta = y(2,:);

%---Beam-Warming

i = 1:length(x);

t = 0;

T = 4/a + dt;

x0 = -6;


A = zeros(length(x),length(x));
    
    id = ones(1,length(x));
    L_zeta = lambda_1*dt/(4*h)*ones(1,length(x)-1);
    R_zeta = -lambda_1*dt/(4*h)*ones(1,length(x)-1);
    
    A_zeta = A + diag(id) + diag(L_zeta,-1) + diag(R_zeta,1);
    
    A_zeta(1,end) = lambda_1*dt/(4*h);
    A_zeta(end,1) = -lambda_1*dt/(4*h);
    
    
    L_eta = lambda_2*dt/(4*h)*ones(1,length(x)-1);
    R_eta = -lambda_2*dt/(4*h)*ones(1,length(x)-1);
    
    A_eta = A + diag(id) + diag(L_eta,-1) + diag(R_eta,1);
    
    A_eta(1,end) = lambda_2*dt/(4*h);
    A_eta(end,1) = -lambda_2*dt/(4*h);
    
    
    [Lz,Uz] = lu(A_zeta);
    [Le,Ue] = lu(A_eta);

while t < T
    
    R = mod(i,length(x))+1;
    L = mod(i-2,length(x))+1;
    
    b_zeta = zeta - 0.5*(dt/h)*(-lambda_1*zeta(R)+lambda_1*zeta(L)) + (dt/(4*h))*(-lambda_1*zeta(R)+lambda_1*zeta(L));
    
    b_eta = eta - 0.5*(dt/h)*(-lambda_2*eta(R)+lambda_2*eta(L)) + (dt/(4*h))*(-lambda_2*eta(R)+lambda_2*eta(L));
    
    b_zeta = b_zeta';
    b_eta = b_eta';
    
    zeta = Uz\(Lz\b_zeta);
    eta = Ue\(Le\b_eta);
    
    zeta = zeta';
    eta = eta';
    
    
    y = [zeta;eta];
    w = C*y;
    
    j = 2:length(x);
    
    u(1) = 0;
    
    sum = 0;
       
    for j = 2:length(x)
        
        sum = sum + 0.5*h*(w(2,j-1) + w(2,j));
        
        u(j) = sum;
        
    end
    
    
    format long
    
    if (abs(t- 1/(2*a))<=1e-12 || abs(t-1/a)<=1e-12 || abs(t-2/a)<=1e-12 || abs(t-4/a)<=1e-12)
        
        figure;
        hold on
        str = sprintf('Beam-Warming, t = %.1f(1/a)',round(t*a,1));
        title(str);
        plot(x,u);
        u_exact = 0.5*(phi_func(x + a*t) + phi_func(x - a*t));
        plot(x,u_exact);
        
        
        
    end
    
    t = t + dt;
    
end







function out = phi_func(x)

    zero = zeros(1,length(x));
    comp = [1 - abs(x);zero];
    out = max(comp,[],1);

end
