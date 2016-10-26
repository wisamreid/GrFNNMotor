function pdex4

m=0; %slab
r=linspace(0,1,100);
t=linspace(1500,1501,3);
sol=pdepe(m,@pdex4pde,@pdex4ic,@bc2fn,r,t);
disp(sol);

u1 = sol(:,:,1);
u2 = sol(:,:,2);

figure
% surf(r,t,u1)
plot(r,u1)
title('X vs r')
xlabel('Distance r')
ylabel('X')

figure
plot(r,u2)
% surf(r,t,u2)
% title('Y(r,t)')
% xlabel('Distance r')
% ylabel('Time t')
title('Y vs r')
xlabel('Distance r')
ylabel('Y')