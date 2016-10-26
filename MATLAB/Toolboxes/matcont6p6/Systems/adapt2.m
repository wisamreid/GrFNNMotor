function out = adapt2
out{1} = @init;
out{2} = @fun_eval;
out{3} = @jacobian;
out{4} = @jacobianp;
out{5} = @hessians;
out{6} = @hessiansp;
out{7} = @der3;
out{8} = [];
out{9} = [];
out{10}= @x1;
out{11}= @x2;
out{12}= @x3;
out{13}= @x4;
out{14}= @x5;

% --------------------------------------------------------------------------
function dydt = fun_eval(t,kmrgd,alpha,beta)
dydt=[kmrgd(2);
kmrgd(3);
-alpha*kmrgd(3)-beta*kmrgd(2)-kmrgd(1)+kmrgd(1)^2;];

% --------------------------------------------------------------------------
function [tspan,y0,options] = init
handles = feval(adapt2);
y0=[0,0,0];
options = odeset('Jacobian',handles(3),'JacobianP',handles(4),'Hessians',handles(5),'HessiansP',handles(6));
tspan = [0 10];

% --------------------------------------------------------------------------
function jac = jacobian(t,kmrgd,alpha,beta)
jac=[ 0 , 1 , 0 ; 0 , 0 , 1 ; 2*kmrgd(1) - 1 , -beta , -alpha ];
% --------------------------------------------------------------------------
function jacp = jacobianp(t,kmrgd,alpha,beta)
jacp=[ 0 , 0 ; 0 , 0 ; -kmrgd(3) , -kmrgd(2) ];
% --------------------------------------------------------------------------
function hess = hessians(t,kmrgd,alpha,beta)
hess1=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 2 , 0 , 0 ];
hess2=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
hess3=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
hess(:,:,1) =hess1;
hess(:,:,2) =hess2;
hess(:,:,3) =hess3;
% --------------------------------------------------------------------------
function hessp = hessiansp(t,kmrgd,alpha,beta)
hessp1=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , -1 ];
hessp2=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , -1 , 0 ];
hessp(:,:,1) =hessp1;
hessp(:,:,2) =hessp2;
%---------------------------------------------------------------------------
function tens3  = der3(t,kmrgd,alpha,beta)
tens31=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens32=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens33=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens34=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens35=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens36=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens37=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens38=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens39=[ 0 , 0 , 0 ; 0 , 0 , 0 ; 0 , 0 , 0 ];
tens3(:,:,1,1) =tens31;
tens3(:,:,1,2) =tens32;
tens3(:,:,1,3) =tens33;
tens3(:,:,2,1) =tens34;
tens3(:,:,2,2) =tens35;
tens3(:,:,2,3) =tens36;
tens3(:,:,3,1) =tens37;
tens3(:,:,3,2) =tens38;
tens3(:,:,3,3) =tens39;
%---------------------------------------------------------------------------
function tens4  = der4(t,kmrgd,alpha,beta)
%---------------------------------------------------------------------------
function tens5  = der5(t,kmrgd,alpha,beta)
function userfun1=x1(t,kmrgd,alpha,beta)
	userfun1=alpha-3
function userfun2=x2(t,kmrgd,alpha,beta)
	userfun2=alpha-4
function userfun3=x3(t,kmrgd,alpha,beta)
	userfun3=0;
function userfun4=x4(t,kmrgd,alpha,beta)
	userfun4=0;
function userfun5=x5(t,kmrgd,alpha,beta)
	userfun5=alpha-6
