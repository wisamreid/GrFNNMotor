function [C,X,Y]=compute_basins()

% This function computes the basins of attraction
% of the 4 fixed points of the syste 

% dx/dt = x - y^3 + 3*sin(x)
% dx/dt = y + x^4 - y^5     
% 
% by using the Newton method

[z1,iter,res] = Newton_sys(@(x) Ffun(x), @(x) Jfun(x),[-0.5 -1]',1e-10,1000);
[z2,iter,res] = Newton_sys(@(x) Ffun(x), @(x) Jfun(x),[-0.2 -0.2]',1e-10,1000);
[z3,iter,res] = Newton_sys(@(x) Ffun(x), @(x) Jfun(x),[1.5 1.5]',1e-10,1000);
[z4,iter,res] = Newton_sys(@(x) Ffun(x), @(x) Jfun(x),[0.1 1]',1e-10,1000);

N=700;

x=linspace(-2,2,N);

[X,Y]=meshgrid(x,x);

C=zeros(size(X));
tol=1e-5;

for i=1:N
    for j=1:N
        [z,iter,res] = Newton_sys(@(x) Ffun(x), @(x) Jfun(x),[X(i,j) Y(i,j)]',1e-10,1000);
        
        if norm(z-z1)<tol
            C(i,j)=1;
        elseif norm(z-z2)<tol
            C(i,j)=2;
        elseif norm(z-z3)<tol  
            C(i,j)=3;
        elseif norm(z-z4)<tol 
            C(i,j)=4;
        end
            
    end
end

fontsize=18;
figure(1)
clf
pcolor(X,Y,C);
colormap jet
shading flat
xlabel('$x_0$','Interpreter','Latex','Fontsize',fontsize)
ylabel('$y_0$','Interpreter','Latex','Fontsize',fontsize)
set(gca,'Fontsize',fontsize,'Xtick',[-2 -1 0 1 2],'Ytick',[-2 -1 0 1 2])
axis([-2 2 -2 2])
            
            
