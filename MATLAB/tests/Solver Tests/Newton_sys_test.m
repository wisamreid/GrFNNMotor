function [z1,z2,z3,z4] = Newton_sys_test()

% This function computes the four zeros of the nonlinear system 
%
% F1(x,y) = x - y^3 + 3*sin(x) = 0
% F2(x,y) = y + x^4 - y^5      = 0
% 
% by using the Newton method for systems of nonlinear equations



    [X Y] = meshgrid(-2:0.1:2,-2:0.1:2);
    Z1 = X - Y.^3+3.0*sin(X);
    Z2 = Y + X.^4 - Y.^5;
    
    figure(1);
    clf;
    hold on
    v=[0,1.5,3];
   
    [C1,h1]= contour(X,Y,Z1,v,'-r');
    set(h1,'ShowText','on');
    
    [C2,h2]=contour(X,Y,Z2,v,'--b');
    set(h2,'ShowText','on');    
    grid on;
    axis equal;
    axis tight;
    xlabel('x');
    ylabel('y');
    
    tol = 1.0e-10;
    nmax = 100;
    x0 = [-0.5;-1.5];
    [z1,F,niter] = Newton_sys(@Ffun,@Jfun,x0,tol,nmax);
    
    x0 = [0.0;0.001];
    [z2,F,niter] = Newton_sys(@Ffun,@Jfun,x0,tol,nmax);
    
    x0 = [0.2;0.75];
    [z3,F,niter] = Newton_sys(@Ffun,@Jfun,x0,tol,nmax);
    
    x0 = [1.5;1.5];
    [z4,F,niter] = Newton_sys(@Ffun,@Jfun,x0,tol,nmax); 
    
    plot(z1(1),z1(2),'ro');
    plot(z2(1),z2(2),'ro');
    plot(z3(1),z3(2),'ro');
    plot(z4(1),z4(2),'ro');
    
    
    return;
end