H = 2;

W = 2*pi;

Nx = 100;
Ny = 200;

dx = W/Nx;
dy = H/Ny;

x = -pi:dx:(pi-dx);

y = dy:dy:H;


diag_id = repelem([-(2/(dx^2)) - (2/(dy^2))],Nx);
diag_E = (1./(dx^2))*repelem([1],Nx-1);
diag_w = (1./(dx^2))*repelem([1],Nx-1);
diag_N = (1./(dy^2))*repelem([1],Nx*(Ny-1));
diag_S = (1./(dy^2))*repelem([1],Nx*(Ny-1));


diag_N(1:Nx) = 2;

A = zeros(Nx) + diag(diag_id,0) + diag(diag_E,1) + diag(diag_w,-1);

A(1,end) = 1;
A(end,1) = 1;


Ar = repmat(A, 1, Ny);
Ac = mat2cell(Ar, size(A,1), repmat(size(A,2),1,Ny));
L = blkdiag(Ac{:});

B = L + diag(diag_N,Nx) + diag(diag_S,-Nx);



%--------------- source ------------------------------


b = zeros(Nx*Ny,1);

for k = 1:length(b)

   i = mod(k-1,Nx)+1; 
   
   x_i = x(i);
   
   if abs(x_i) <= pi/2
       
       b(k) = -cos(x_i);
       
   end
    
    
end


disp(b);

%------------------------- solver ------------------------------


[L,U] = lu(B);

v = L\b;

u = U\v;

u = reshape(u,Nx,Ny);

u = u';

disp(size(u));


[X,Y] = meshgrid(x,y);

disp(size(X));

surf(X, Y, u);
shading interp
view(0,90);
colormap(jet);
xlabel("x")
ylabel("y")
title("Stationary Heat Distribution (u)");