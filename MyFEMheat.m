function MyFEMheat()
%close all
c = imread('cat.png');
cc = sum(c,3);
h = contour(cc,[1 1]);
% extract contours of face and eyes
ind = find(h(1,:)==1);
c1 = h(:,2:6:ind(2)-1)';  % face
xmin = min(c1(:,1));
xmax = max(c1(:,1));
ymin = min(c1(:,2));
ymax = max(c1(:,2));
midpt = 0.5*[(xmin+xmax),(ymin+ymax)];
dd = 0.5*[(-xmin+xmax),(-ymin+ymax)];
c2 = h(:,ind(3)+1:6:ind(4)-1)'; % eye
c3 = h(:,ind(4)+1:6:ind(5)-1)'; % another eye
lc1 = length(c1); 
lc2 = length(c2);
lc3 = length(c3);
% connect boundary points
a1 = [1 : lc1]';
e1 = [a1, circshift(a1,[-1 0])];
a2 = [1 : lc2]';
e2 = [lc1 + a2, lc1 + circshift(a2,[-1 0])];
a3 = [1 : lc3]';
e3 = [lc1 + lc2 + a3, lc1 + lc2 + circshift(a3,[-1 0])];
c1 = rescale2unitsquare(c1,midpt,dd);
c2 = rescale2unitsquare(c2,midpt,dd);
c3 = rescale2unitsquare(c3,midpt,dd);
bdry = [c1; c2; c3]; % coordinates of boundary points
bdry_connect = [e1; e2; e3]; % indices of endpoints of boundary segments
% load mesh generated using D. Engwirda's routine
msh = load('MyFEMcat_mesh.mat');
pts = msh.pts;
tri = msh.tri;
pts = rescale2unitsquare(pts,midpt,dd);

%---------------------------------------------- do mesh-opt.

% pts is a N-by-2 array with coordinated of the mesh points
% tri is a Ntriag-by-3 array of indices of the triangular elements


% Find boundary points with homogeneous Neumann BCs
% The number of rows in neumann is the number of boundary intervals with
% Neumann BCs

% Each row of neumann contains indices of endpoints of the corresponding
% boundary interval
ind = zeros(lc1,1);
for i = 1 : lc1
    ind(i) = find(pts(:,1) == c1(i,1) & pts(:,2) == c1(i,2));
end
neumann = [ind,circshift(ind,[-1, 0])];

% Find boundary points with Dirichlet BCs
% dirichlet is a column vector of indices of points with Dirichlet BCs
dirichlet = zeros(lc2 + lc3,1);
ind = zeros(lc2,1);
for i = 1 : lc2
    dirichlet(i) = find(pts(:,1) == c2(i,1) & pts(:,2) == c2(i,2));
end
ind = zeros(lc3,1);
for i = 1 : lc3
    dirichlet(lc2 + i) = find(pts(:,1) == c3(i,1) & pts(:,2) == c3(i,2));
end

% call FEM
fem2d_heat(pts,tri,neumann,dirichlet);
end



function fem2d_heat(pts,tri,neumann,dirichlet)
%FEM2D_HEAT finite element method for two-dimensional heat equation. 2 %Initialisation
FreeNodes=setdiff(1:size(pts,1),unique(dirichlet));
Npts = size(pts,1);
A = sparse(Npts,Npts);
B = sparse(Npts,Npts); 
T = 1; 
dt = 0.01; 
N = T/dt;
tdraw = 0.2;
jdraw = tdraw/dt;
Ndraw = T/tdraw+1;
% Assembly
for j = 1:size(tri,1)
	A(tri(j,:),tri(j,:)) = A(tri(j,:), ...
    tri(j,:)) + stima3(pts(tri(j,:),:));
end
for j = 1:size(tri,1)
	B(tri(j,:),tri(j,:)) = B(tri(j,:), ...
    tri(j,:)) + det([1,1,1;pts(tri(j,:),:)'])... 
    *[2,1,1;1,2,1;1,1,2]/24;
end
%Initial Condition
u = IC(pts); 
plotsolution(pts,tri,u,0)
%save the initial temperature distribution
Udraw = zeros(Npts,Ndraw);
Udraw(:,1) = u;
%Set up the movie.
fps = 24; %frames per second -- normally 24 frames per second
writerObj = VideoWriter('MyFEMheat.mp4','MPEG-4'); % Name it.
writerObj.FrameRate = fps; % How many frames per second.
open(writerObj); 
iframe = 0;
% get the first frame
frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
writeVideo(writerObj, frame);
iframe = iframe + 1;
drawcount = 1;
 
% time steps
for n = 2:N+1
    uold = u;
    b = sparse(size(pts,1),1);
    t = (n-1)*dt; % time at the new step
    % Volume Forces
    for j = 1:size(tri,1)
        b(tri(j,:)) = b(tri(j,:)) + ... 
            det([1,1,1; pts(tri(j,:),:)']) * ... 
            dt*myf(sum(pts(tri(j,:),:))/3,t+dt/2)/6;
    end
    % Neumann conditions
    for j = 1 : size(neumann,1)
       b(neumann(j,:)) = b(neumann(j,:)) + ...
         norm(pts(neumann(j,1),:)-pts(neumann(j,2),:))*...
         dt*myg(sum(pts(neumann(j,:),:))/2,t+dt/2)/2;
    end
    b=b+B*uold-(dt/2)*A*uold;
    % Homogeneous Dirichlet conditions
    u = sparse(size(pts,1),1);
    u(unique(dirichlet)) = zeros(size(unique(dirichlet))); 
    % Computation of the solution
    u(FreeNodes) = ((dt/2)*A(FreeNodes,FreeNodes)+ ...
            B(FreeNodes,FreeNodes))\b(FreeNodes);
    plotsolution(pts,tri,u,t);
    frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
    writeVideo(writerObj, frame);
    iframe = iframe + 1;
    % save u for future drawing
    if mod(n,jdraw) == 0
        drawcount = drawcount+1;
        U(:,drawcount) = u;
    end
end % end of time marching

close(writerObj); % Saves the movie.
%plot figures
for k = 1 : drawcount
    figure;
    t = 0.2*(k - 1);
    u = U(:,k);
    plotsolution(pts,tri,u,t);
end
end

function u0 = IC(x)
u0 = ones(size(x,1),1); 
end


function Stress = myg(x,~)
Stress = zeros(size(x,1),1);
end


function M = stima3(vertices)
d = size(vertices,2);
G = [ones(1,d+1);vertices'] \ [zeros(1,d);eye(d)];
M = det([ones(1,d+1);vertices']) * G * G' / prod(1:d);
end


function heatsource = myf(x,t)
heatsource = 10*exp(-4*(x(:,1).^2+(x(:,2)-0.8).^2))*cos(pi*t);
end

function x = rescale2unitsquare(x,midpt,dd)
% x should be an array n-by-2
e = ones(size(x,1),1);
x = (x-e*midpt)./(e*dd);
end

function plotsolution(pts,tri,u,t)
trisurf(tri,pts(:,1),pts(:,2),full(u)','facecolor','interp')
axis ij
title(sprintf('Time = %.4f',t),'Fontsize',24);
colorbar;
caxis([-1,1]);
set(gca,'Fontsize',20);
view(2)
end

