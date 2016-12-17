function compute_bifurcation_diagram_1D(f,df,z0,z1,x0,x1)
% This function computes the bifurcation diagram (equilibra)  
% for an ODE in the form 
%
% dx/dt = f(x;z) - z is the bifurcation parameter
%
% Input  f -> function handle in the form f=@(x,z) x^2*z+sin(z*x)
%       df -> derivative f(x,z) with respect to x 
%       [z0,,z1] -> lower and upper bound on the bifurcation parameter
%       [x0,x1]  -> lower and upper bound on the fixed points we are i
%                   interested in

Nz=400;  % This is the number of points between z0 and z1
Nx= 15; % This is the number of initial guesses for our solver


Z=linspace(z0,z1,Nz);
X=linspace(x0,x1,Nx);


fontsize=14;
figure(1)
clf
hold
box on
color =['b' 'r' 'k'];
for ii=1:Nz
    
    for jj=1:Nx
       
        %[zz,fval,flag]=fzero(@ (x,z) f(x,Z(ii)), X(jj));
  [zz,iter,res] = Newton_method(@ (x,z) f(x,Z(ii)),@ (x,z) df(x,Z(ii)),X(jj),1e-5,100);
        if res < 1e-4
            idx=check_stability(f,zz,Z(ii));
          %  linetype=sprintf('%s.',color(idx));
          %  plot(Z(ii),zz,linetype);
          scatter(Z(ii),zz,4,idx,'filled')
        end
    end
end

axis([z0 z1 x0 x1])
xlabel('$\mu$', 'Interpreter','Latex','Fontsize',fontsize);
ylabel('$x$', 'Interpreter','Latex','Fontsize',fontsize);
set(gca,'Fontsize',fontsize);
grid


function idx=check_stability(f,x,z)

% Output idx -> index that tells us if the node is stable (1), unstable (2) 
%               or half-stable (3)
tol=1e-2;
a1=f(x+tol,z); 
a2=f(x-tol,z);

if  (a1*a2 <0) && (a2>0) % Stable node
    idx=1;
elseif (a1*a2<0) && (a2< 0) % Unstable node
    idx=2;
elseif (a1*a2 > 0)  % Half stable 
    idx=3;
end