

[A,b] = MyFEMp1();

L = ichol(A);

L_inv = L\eye(size(L));

A_t = L_inv*A*(L_inv)';

b_t = L_inv*b;


u_0 = zeros(length(A_t),1);

u = u_0;

r = b_t - A*u_0;

p = r;

res = r;

res_list = [];

res_list(end+1) = vecnorm(res);

while vecnorm(res) > 1e-12
    
    alpha = (r')*r/((p')*A_t*p);
    
    u = u + alpha*p;
    r_new = r - alpha*A_t*p;
    
    beta = (r_new')*r_new/((r')*r);
    
    p = r_new + beta*p;
    
    r = r_new;
    
    res = r;
    
    res_list(end+1) = vecnorm(res);
    
    
end 

x = (L_inv')*u;

figure;
semilogy(0:(length(res_list)-1),res_list);
grid on;
title('log(residual norm)');
xlabel('k (iteration number)');
ylabel('log(||b - Au||)');

