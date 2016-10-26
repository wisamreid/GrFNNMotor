function [y,ts] = AB3_method(F,y0,DT,N,IOSTEPS)

% This function computes the solution to dy/dt=F(y,t) y(t0)=y0 in the 
% interval [0,T] by using the Adams-Bashforth 3th order method

% Input    F  -> Function handle characterizing F(y,t)
%          y0 -> Initial condition vector (Has to be a column vector)
%          DT -> Delta t
%           N -> Total number of steps 
%     IOSTEPS -> We dump out a snapshot every IOSTEPS steps




SNAPS   = floor(N/IOSTEPS)+1;

y=zeros(size(y0,1),SNAPS);

T1=N*DT;

t=linspace(0,DT*N,N+1);


y(:,1)=y0; % This is the initial condition

% We compute y1 and y2 required to start-up the Adams Bashforth method
% by using the Heun method

y1 = y0 + DT/2* ( F(y0+DT*F(y0,t(1)),t(2)) +F(y0,t(1)));
y2 = y1 + DT/2* ( F(y1+DT*F(y1,t(2)),t(3)) +F(y1,t(2)));

I=1;

for i=3:N

% Hereafter is RK4 explicit scheme
f2 = F(y2,t(i));
f1 = F(y1,t(i-1));
f0 = F(y0,t(i-2));

y3= y2 + DT/12* (23*f2-16*f1+5*f0);

if mod(i,IOSTEPS)==0
   I=I+1;
   ts(I)=t(i);
   y(:,I)=y3;
end

y0=y1;
y1=y2;
y2=y3;

end

end




