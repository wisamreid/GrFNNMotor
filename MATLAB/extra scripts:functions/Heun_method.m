function [y,ts] = Heun_method(F,y0,DT,N,IOSTEPS)

% This function computes the solution to dy/dt=F(y,t) y(t0)=y0 in the 
% interval [0,T] by using the Heun method

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

I=1;
for i=1:N

y1= y0 + DT/2* ( F(y0+DT*F(y0,t(i)),t(i+1)) +F(y0,t(i)));

if mod(i,IOSTEPS)==0
   I=I+1;
   ts(I)=t(i);
   y(:,I)=y1;
end

y0=y1;

end

end

