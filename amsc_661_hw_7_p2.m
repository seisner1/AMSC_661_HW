
%-------------------- problem 2
%---- code link: 

Ms = [1,2,5,10,100];
X = (0:0.01:pi);
close all;

for i=1:5
    
    m = Ms(i);
    
    U = u(m);
    
    figure;
    plot(X,U(X,0));
    xlabel('x');
    ylabel('u(x,0)');
    str = sprintf('sum of first %d terms of u(x,0)',m);
    title(str);
    
    figure;
    plot(X,U(X,2));
    xlabel('x');
    ylabel('u(x,2)');
    str = sprintf('sum of first %d terms of u(x,2)',m);
    title(str);
    
    
    if m==100
        
        max_norm = max(abs(X - U(X,0)));
        
        max_norm
        
    end
    
    
    
end

function out = u(m)

    N = (0:1:m)';
    
    out = @(x,t) sum((2.*((-1).^N)./(pi.*(N+0.5).^2)).*sin((N+0.5)*x).*exp(-t.*(N+0.5).^2),1);

end