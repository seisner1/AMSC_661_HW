
function [A,b] = MyFEMp1()
close all
fd = inline('dunion(dunion(dcircle(p,-1,3,4),dcircle(p,-3,3,0.5)),dcircle(p,0,4.5,0.5))','p');

p_rect = [0,0;3,0;3,3;0,3];
phi_circ = (0:99)'/100*2*pi + pi/2;

phi_circ_d = (0:19)'/20*2*pi + pi/2;

p_gamma = 4.*[cos(phi_circ),sin(phi_circ)] + [-1,3];

p_a = 0.5.*[cos(phi_circ_d),sin(phi_circ_d)] + [-3,3];

p_b = 0.5.*[cos(phi_circ_d),sin(phi_circ_d)] + [0,4.5];

p_fix = [p_gamma;p_a;p_b];

[pts,tri] = distmesh2d(fd,@huniform,0.05,[-5,3;-1,7],p_fix);

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
ind = find(abs((pts(:,1) + 1).^2 + (pts(:,2) - 3).^2 - 16) < 0.05);
ind = ind;

neumann = [ind,circshift(ind,[-1, 0])];

% Find boundary points with Dirichlet BCs
% dirichlet is a column vector of indices of points with Dirichlet BCs
dirichlet = zeros(lc1 + lc1,1);
ind = zeros(lc1,1);
dirichlet = find(abs((pts(:,1) + 3).^2 + (pts(:,2) - 3).^2 - 0.25) < 0.05);
dirichlet = [dirichlet;find(abs((pts(:,1) - 0).^2 + (pts(:,2) - 4.5).^2 - 0.25) < 0.05)];

% call FEM
[u,A_o,b_o] = FEM2D(pts,tri,neumann,dirichlet,1.2,1);

% graphic representation

figure;
trisurf(tri,pts(:,1),pts(:,2),full(u),'facecolor','interp')
colormap(jet);
hold on
axis ij
view(2)
xlabel('x','Fontsize',20);
ylabel('y','Fontsize',20);
set(gca,'Fontsize',20, 'YDir', 'normal')
colorbar
title('u(x,y)');

A = A_o;
b = b_o;

% figure;
% 
% current_dens = zeros(size(tri,1),1);
% for i = 1:size(tri,1)
%     
%     verts = pts(tri(i,:),:);
%     
%     mdpt = [(1/3)*sum(verts(:,1)),(1/3)*sum(verts(:,2))];
%     
%     G = [ones(1,3);verts'] \ [zeros(1,2);eye(2)];
%    
%     grad_u = u(tri(i,:))'*G;
%     
%     
%     current_dens(i) = norm(-a(mdpt,1.2,1)*grad_u);
%     
% end
% 
% Npts = size(pts,1);
% Ntri = size(tri,1);
% 
% abs_current_verts = zeros(Npts,1);
% count_tri = zeros(Npts,1);
% for j = 1:size(tri,1)
% abs_current_verts(tri(j,:)) = abs_current_verts(tri(j,:)) ...
% + current_dens(j);
% count_tri(tri(j,:)) = count_tri(tri(j,:)) + 1;
% end
% abs_current_verts = abs_current_verts./count_tri;
% 
% trisurf(tri,pts(:,1),pts(:,2),full(abs_current_verts),'facecolor','interp')
% colormap(jet);
% hold on
% axis ij
% view(2)
% xlabel('x','Fontsize',20);
% ylabel('y','Fontsize',20);
% set(gca,'Fontsize',20);
% colorbar
% title('Abs. Current Density, a_1=1.2, a_2 = 1');


% %------------------------------a_1 = 0.8---------------
% 
% % call FEM
% u = FEM2D(pts,tri,neumann,dirichlet,1.2,1);
% 
% % graphic representation
% figure;
% trisurf(tri,pts(:,1),pts(:,2),full(u),'facecolor','interp')
% colormap(jet);
% hold on
% axis ij
% view(2)
% xlabel('x','Fontsize',20);
% ylabel('y','Fontsize',20);
% set(gca,'Fontsize',20, 'Ydir', 'Reverse');
% colorbar
% title('u(x,y)');

% figure;
% 
% current_dens = zeros(size(tri,1),1);
% for i = 1:size(tri,1)
%     
%     verts = pts(tri(i,:),:);
%     
%     mdpt = [(1/3)*sum(verts(:,1)),(1/3)*sum(verts(:,2))];
%     
%     G = [ones(1,3);verts'] \ [zeros(1,2);eye(2)];
%    
%     grad_u = u(tri(i,:))'*G;
%     
%     
%     current_dens(i) = norm(-a(mdpt,0.8,1)*grad_u);
%     
% end
% 
% Npts = size(pts,1);
% Ntri = size(tri,1);
% 
% abs_current_verts = zeros(Npts,1);
% count_tri = zeros(Npts,1);
% for j = 1:size(tri,1)
% abs_current_verts(tri(j,:)) = abs_current_verts(tri(j,:)) ...
% + current_dens(j);
% count_tri(tri(j,:)) = count_tri(tri(j,:)) + 1;
% end
% abs_current_verts = abs_current_verts./count_tri;
% 
% trisurf(tri,pts(:,1),pts(:,2),full(abs_current_verts),'facecolor','interp')
% colormap(jet);
% hold on
% axis ij
% view(2)
% xlabel('x','Fontsize',20);
% ylabel('y','Fontsize',20);
% set(gca,'Fontsize',20);
% colorbar
% title('Abs. Current Density, a_1=0.8, a_2 = 1');

end

% FEM
function [u,A_o,b_o] = FEM2D(pts,tri,neumann,dirichlet,a_1,a_2)
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

b_o = b;
A_o = A;

end

%
function DirichletBoundaryValue = myu_d(x)
d = zeros(size(x,1),1);
d((round(size(x/1)/2)+1):end) = 1;
DirichletBoundaryValue =  d;
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
x = pt(1);
y = pt(2);
xa=-3; ya=3;
xb=0; yb=4.5;
fac=10;
f=(1-x).^2+(y-0.25*x.^2).^2+1;
g1=1-exp(-0.125*((x-xa).^2+(y-ya).^2));
g2=1-exp(-0.25*(((x-xb).^2+(y-yb).^2)));
g3=1.2-exp(-2.*((x+0).^2+(y-2).^2));
g4=1+exp(-2*(x+1.5).^2-(y-3.5).^2-(x+1).*(y-3.5));
v1=f.*g1.*g2.*g3.*g4;
V=fac*atan(v1/fac);

a = exp(-V);

end