function  [ts,QS]=solver_Lorenz_MC()
% This function solves the Lorenz 96 using ode45.
% The initial state is random 
% (jointly Gaussian with zero mean small std) 
% 
% Output 
%        ts -> vector of time instants 
%        QS -> 3D matrix QS(i,j,k)  solution sample j of xk(ti) k=1,2,3 
%
global R S B

st    = 1e-5;  % (standard deviation of the initial state)

R    = 28;   % Rayleigh number 
S    = 10;   % Prandtl number 
B    = 8/3;

NSAMP = 1;   % number of solution samples

SNAPS =15000;  % Number of time snapshots

tt = linspace(0,200,SNAPS+1); % Time mesh

y0=randn(3,NSAMP);

fprintf('\n ')

for ii=1:NSAMP

tic;

Y0=st*y0(:,ii);  % Initial condition

% options = odeset('RelTol', 1e-4,'AbsTol',1e-4);
[ts,sol]=ode45(@(t,y) RHS(t,y),tt,Y0);

QS(:,ii,:)=sol;
ttt=toc;
fprintf('\n Sample (%d) of (%d)   CPU-TIME=%3.3f',ii,NSAMP,ttt)

end

eval(sprintf('save MC_data_R%3.3g_S%3.3g_B%3.3g_NSAMP%d.mat ts QS ts',R,S,B,NSAMP))


% SAMPLE FIGURE
figure(1)
clf
plot3(QS(:,1,1),QS(:,1,2),QS(:,1,3))
grid
axis tight
xlabel('$x$','Interpreter','Latex','Fontsize',20)
ylabel('$y$','Interpreter','Latex','Fontsize',20)
zlabel('$z$','Interpreter','Latex','Fontsize',20)
set(gca,'Fontsize',20)
view(20,30)

fprintf('\n\n')


function F=RHS(t,x)

% This is the RHS of the Lorenz model
global S R B

F=zeros(3,1);

F(1) = -S*x(1) + S*x(2);
F(2) = R*x(1)-x(2)-x(1)*x(3);
F(3) = x(1)*x(2)-B*x(3);






