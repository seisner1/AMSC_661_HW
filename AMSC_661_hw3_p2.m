function [A,b] = MyFEMCat()
close all
fd = inline('dunion(drectangle(p,0,3,0,3),dcircle(p,1.5,1.5,1))','p');

p_rect = [0,0;3,0;3,3;0,3];
phi_circ = (0:99)'/100*2*pi + pi/2;

p_circ = [cos(phi_circ),sin(phi_circ)] + [1.5,1.5];

p_fix = [p_rect;p_circ];

[pts,tri] = distmesh2d(fd,@huniform,0.1,[0,0;3,3],p_fix);

lc1 = 3/0.1;
%---------------------------------------------- do mesh-opt.

% pts is a N-by-2 array with coordinated of the mesh points
% tri is a Ntriag-by-3 array of indices of the triangular elements

%
% Find boundary points with homogeneous Neumann BCs
% The number of rows in neumann is the number of boundary intervals with
% Neumann BCs

% Each row of neumann contains indices of endpoints of the corresponding
% boundary interval
ind = zeros(lc1 + lc1,1);
ind = find(abs(pts(:,2)) < 0.05);
ind = [ind;find(abs(pts(:,2) - 3) < 0.05)];

neumann = [ind,circshift(ind,[-1, 0])];

% Find boundary points with Dirichlet BCs
% dirichlet is a column vector of indices of points with Dirichlet BCs
dirichlet = zeros(lc1 + lc1,1);
ind = zeros(lc1,1);
dirichlet = find(abs(pts(:,1)) < 0.05);
dirichlet = [dirichlet;find(abs(pts(:,1) - 3) <0.05)];

% call FEM
u = FEM2D(pts,tri,neumann,dirichlet,1.2,1);

% graphic representation
figure;
trisurf(tri,pts(:,1),pts(:,2),full(u),'facecolor','interp')
colormap(jet);
hold on
axis ij
view(2)
xlabel('x','Fontsize',20);
ylabel('y','Fontsize',20);
set(gca,'Fontsize',20);
colorbar
title('Voltage, a_1=1.2, a_2 = 1');

figure;

current_dens = zeros(size(tri,1),1);
for i = 1:size(tri,1)
    
    verts = pts(tri(i,:),:);
    
    mdpt = [(1/3)*sum(verts(:,1)),(1/3)*sum(verts(:,2))];
    
    G = [ones(1,3);verts'] \ [zeros(1,2);eye(2)];
   
    grad_u = u(tri(i,:))'*G;
    
    
    current_dens(i) = norm(-a(mdpt,1.2,1)*grad_u);
    
end

Npts = size(pts,1);
Ntri = size(tri,1);

abs_current_verts = zeros(Npts,1);
count_tri = zeros(Npts,1);
for j = 1:size(tri,1)
abs_current_verts(tri(j,:)) = abs_current_verts(tri(j,:)) ...
+ current_dens(j);
count_tri(tri(j,:)) = count_tri(tri(j,:)) + 1;
end
abs_current_verts = abs_current_verts./count_tri;

trisurf(tri,pts(:,1),pts(:,2),full(abs_current_verts),'facecolor','interp')
colormap(jet);
hold on
axis ij
view(2)
xlabel('x','Fontsize',20);
ylabel('y','Fontsize',20);
set(gca,'Fontsize',20);
colorbar
title('Abs. Current Density, a_1=1.2, a_2 = 1');


%------------------------------a_1 = 0.8---------------

% call FEM
u = FEM2D(pts,tri,neumann,dirichlet,1.2,1);

% graphic representation
figure;
trisurf(tri,pts(:,1),pts(:,2),full(u),'facecolor','interp')
colormap(jet);
hold on
axis ij
view(2)
xlabel('x','Fontsize',20);
ylabel('y','Fontsize',20);
set(gca,'Fontsize',20);
colorbar
title('Voltage, a_1=0.8, a_2 = 1');

figure;

current_dens = zeros(size(tri,1),1);
for i = 1:size(tri,1)
    
    verts = pts(tri(i,:),:);
    
    mdpt = [(1/3)*sum(verts(:,1)),(1/3)*sum(verts(:,2))];
    
    G = [ones(1,3);verts'] \ [zeros(1,2);eye(2)];
   
    grad_u = u(tri(i,:))'*G;
    
    
    current_dens(i) = norm(-a(mdpt,0.8,1)*grad_u);
    
end

Npts = size(pts,1);
Ntri = size(tri,1);

abs_current_verts = zeros(Npts,1);
count_tri = zeros(Npts,1);
for j = 1:size(tri,1)
abs_current_verts(tri(j,:)) = abs_current_verts(tri(j,:)) ...
+ current_dens(j);
count_tri(tri(j,:)) = count_tri(tri(j,:)) + 1;
end
abs_current_verts = abs_current_verts./count_tri;

trisurf(tri,pts(:,1),pts(:,2),full(abs_current_verts),'facecolor','interp')
colormap(jet);
hold on
axis ij
view(2)
xlabel('x','Fontsize',20);
ylabel('y','Fontsize',20);
set(gca,'Fontsize',20);
colorbar
title('Abs. Current Density, a_1=0.8, a_2 = 1');

end

% FEM
function u = FEM2D(pts,tri,neumann,dirichlet,a_1,a_2)
Npts = size(pts,1);
Ntri = size(tri,1);
FreeNodes = setdiff(1:Npts,dirichlet); %mesh points with unknown values of u
A = sparse(Npts,Npts);
b = sparse(Npts,1);

% Assembly
% The Stiffness matrix
for j = 1:Ntri % for all triangles    
  A(tri(j,:),tri(j,:)) = A(tri(j,:),tri(j,:)) + stima3(pts(tri(j,:),:),a_1,a_2);
  % stima3_2 computes M = 0.5*|T_j|*G*G';
end

% The Right-hand side, i.e., the load vector
% Volume Forces
for j = 1:Ntri
  b(tri(j,:)) = 0;  % for the case where f = 0
end

% Neumann conditions
for j = 1 : size(neumann,1)
  b(neumann(j,:)) = b(neumann(j,:)) + norm(pts(neumann(j,1),:)- ...
      pts(neumann(j,2),:)) * myg(sum(pts(neumann(j,:),:))/2)/2;
end

% Dirichlet conditions 
u = sparse(Npts,1);
u(dirichlet) = myu_d(pts(dirichlet,:));
b = b - A * u;

% Computation of the solution
u(FreeNodes) = A(FreeNodes,FreeNodes) \ b(FreeNodes);

b = b;
A = A;

end

%
function DirichletBoundaryValue = myu_d(x)
xmin = min(x(:,1));
xmax = max(x(:,1));
midx = 0.5*(xmin + xmax);
DirichletBoundaryValue =  0.5 * (sign(x(:,1) - midx) + 1);
end

%
function Stress = myg(x)
Stress = zeros(size(x,1),1);
end

%
function M = stima3(verts,a_1,a_2)
G = [ones(1,3);verts'] \ [zeros(1,2);eye(2)];
mdpt = [(1/3)*sum(verts(:,1)),(1/3)*sum(verts(:,2))];
M = a(mdpt,a_1,a_2)*0.5*det([ones(1,3);verts']) * G * G';
end

function a = a(pt,a_1,a_2)
if ((pt(1)-1.5)^2 + (pt(2)-1.5)^2 <= 1)
    
    a = a_1;
else
    
    a = a_2;
end

end