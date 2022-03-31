% N = [6 11];

%-----------------------problem 1.b------
%----- code link: 

close all
times = [];
Ns = [];
for i = 25:1:100
    N = i.^2;
%     fprintf('\n');
    tstart = cputime;
    NestedDissection(i);
    tend = cputime;
    
    t = tend - tstart;
    
    times(end+1) = t;
    Ns(end+1) = N;
    
end

figure;
grid on
loglog(Ns,times);
title('log of CPU time vs. log of N = n0^2');
bool = (~isnan(log(times))) & (~isnan(log(Ns)));
c = polyfit(log(Ns(6:end)),log(times(6:end)),1);

fprintf('p = %f \n',c(1));

function NestedDissection(n0)
% n0 must be of the form 2^{k+1} - 1 or 2^{k+1} + 2^k - 1
% so that the matrices A11 and A22 are of equal sizes
% n0 = 3, 7, 15, 31, 63, 127 -- do not try larger, or
% n0 = 5, 11, 23, 47, 95
[A,b,sol] = TestMatrixA(n0 + 2);
n = n0;
level = 0;
[L,U,P,A] = MyDissection(A,n,n,level);
y = L\(P*b);
x = P'*(U\y);

% fprintf('norm(x - sol) = %d\n',norm(x - sol));
% figure;
% [I,J] = ind2sub(size(A),find(A));
% plot(I,J,'.','Markersize',10);
% grid
% title(sprintf('Sparsity pattern of P*A*P^T for nx = %d, ny = %d',n0,n0),'Fontsize',20);
% axis ij

end


function [L,U,P,A] = MyDissection(A,nx,ny,level)
A0 = A;
[m,n] = size(A);
if m ~= n
    fprintf("A is not square: (%d-by-%d)\n",m,n);
    L=0;
    U=0;
    P=0;
    A=0;
    return
end
% if level is even do vertical split
% if level is odd do horizontal split
nxy = nx*ny;
par = mod(level,2);% parity
e_flag = 0;
switch par
    case 0 % vertical split
        if nx >= 3
            m = floor(nx/2);
            mn = m*ny;
%             ind = mn + 1 : mn + ny; % indices of Omega3
%             ind1 = 1 : mn; % indices of Omega1
%             ind2 = mn + ny + 1 : nxy; % indices of Omega2
%             nxnext = m;
%             nynext = ny;
%             disp(ind1);
            if mod(nx,2)~= 0
                ind1 = 1 : mn; % indices of Omega1
                ind2 = mn + ny + 1:nxy;
                ind = mn + 1 : mn + ny; % indices of Omega3
                nxnext = m;
                nynext = ny;
                e_flag = 0;
            end
%             
            if mod(nx,2)==0
               ind1 = 1 : (mn - ny); % indices of Omega1
               ind2 = mn + ny + 1:nxy;
               ind = mn - ny + 1 : mn + ny; % indices of Omega3
               nxnext = m-1;
               nynext = ny;
               e_flag = 1;
            end

        else
            [L,U] = lu(A);
            P = speye(nxy);
            return
        end    
    case 1 % horizontal split
        if ny >= 3
            m = floor(ny/2);
            mn = m*nx;
%             ind = m + 1 : ny : nxy; % indices of Omega3
%             [ii,jj] = meshgrid(1 : m,1 : nx);
%             ind1 = sort(sub2ind([ny,nx],ii(:)',jj(:)'),'ascend'); % indices of Omega1
%             [ii,jj] = meshgrid(m + 2 : ny,1 : nx);
%             ind2 = sort(sub2ind([ny,nx],ii(:)',jj(:)'),'ascend'); %indices of Omega2 
%             nxnext = nx;
%             nynext = m;
            
            if mod(ny,2)~= 0
                ind = m + 1 : ny : nxy; % indices of Omega3
                [ii,jj] = meshgrid(1 : m,1 : nx);
                ind1 = sort(sub2ind([ny,nx],ii(:)',jj(:)'),'ascend'); % indices of Omega1
                [ii,jj] = meshgrid(m + 2 : ny,1 : nx);
                ind2 = sort(sub2ind([ny,nx],ii(:)',jj(:)'),'ascend'); %indices of Omega2 
                nxnext = nx;
                nynext = m;
                e_flag = 0;
                
            end
            
            if mod(ny,2)==0
                [ii,jj] = meshgrid(m:(m+1),1:nx);
                ind = sort(sub2ind([ny,nx],ii(:)',jj(:)'),'ascend'); % indices of Omega3
                [ii,jj] = meshgrid(1 : (m-1),1 : nx);
                ind1 = sort(sub2ind([ny,nx],ii(:)',jj(:)'),'ascend'); % indices of Omega1
                [ii,jj] = meshgrid(m + 2 : ny,1 : nx);
                ind2 = sort(sub2ind([ny,nx],ii(:)',jj(:)'),'ascend'); %indices of Omega2 
                nxnext = nx;
                nynext = m-1;
                e_flag = 1;
                
            end
        else
            [L,U] = lu(A);
            P = speye(nxy);
            return
        end    
    otherwise
%         fprintf('Error: par = %d\n',par);
        return
end
A11 = A(ind1,ind1);
A22 = A(ind2,ind2);

[L11,U11,P11,A11] = MyDissection(A11,nxnext,nynext,level + 1);
[L22,U22,P22,A22] = MyDissection(A22,nxnext,nynext,level + 1);

P1 = speye(nxy);
P1(ind1,ind1) = P11;
P1(ind2,ind2) = P22;
% set up the permutation matrix P
P = sparse(1 : nxy,[ind1(:)',ind2(:)',ind(:)'],ones(1,nxy));
P = P*P1;
A = P*A0*P';
% extract nonzero blocks of A
A11 = A(ind1,ind1);
mn1 = mn + 1;
mn2 = 2*mn;
mn3 = mn2 + 1;
A22 = A(mn1:mn2,mn1:mn2);
A13 = A(1 : mn,mn3 : end);
A23 = A(mn1 : mn2,mn3 : end);
A31 = A(mn3 : end,1 : mn);
A32 = A(mn3 : end,mn1 : mn2);
A33 = A(mn3:end,mn3:end);

if e_flag
    
    if (par==0)
    nn = ny;
    end
    if (par ==1)
        nn = nx;
    end
   
    A11 = A(ind1,ind1);
    mn1 = mn - nn + 1;
    mn2 = 2*(mn - nn);
    mn3 = mn2 + 1;
    A22 = A(mn1:mn2,mn1:mn2);
    A13 = A(ind1,mn3 : end);
    A23 = A(mn1 : mn2,mn3 : end);
    A31 = A(mn3 : end,ind1);
    A32 = A(mn3 : end,mn1 : mn2);
    A33 = A(mn3:end,mn3:end);
    
    
    
end
% compute the Schur compliment
S33 = A33 - A31*(U11\(L11\A13)) - A32*(U22\(L22\A23));
% compute LU factorization of S33
[L33,U33] = lu(S33);
% form the LU decomposition of A
L = sparse(nxy,nxy);
first_ind = 1:mn;
if e_flag
    
    if (par==0)
    nn = ny;
    end
    if (par ==1)
        nn = nx;
    end
   
    first_ind = 1:(mn-nn);
    mn1 = mn - nn + 1;
    mn2 = 2*(mn - nn);
    mn3 = mn2 + 1;
    
end
L(first_ind,first_ind) = L11;
L(mn1 : mn2,mn1 : mn2) = L22;
L(mn3 : end,mn3 : end) = L33;
L(mn3 : end,first_ind) = (U11'\A31')';
L(mn3 : end,mn1 : mn2) = (U22'\A32')';
U = sparse(nxy,nxy);
U(first_ind,first_ind) = U11;
U(mn1 : mn2,mn1 : mn2) = U22;
U(mn3 : end,mn3 : end) = U33;
U(first_ind,mn3 : end) = L11\A13;
U(mn1 : mn2,mn3 : end) = L22\A23;

% fprintf('nx = %d, ny = %d, level = %d, size(A) = [%d,%d]\n',nx,ny,level,n,n);
% fprintf('norm(full(A - L*U)) = %d\n',norm(full(A - L*U)));
end

        
        
        
        
        

        
