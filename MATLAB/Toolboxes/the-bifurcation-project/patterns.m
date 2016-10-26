A = 2;
B = 4;
Bc = 3.602;

L = 1;
beta = 631.0144;
Dx = 1.6*10^-3;
nc = 8;
fi = -3/4*(beta*Dx + 1-Bc)/(A^2) + 1/A^2*(2*(beta*Dx + 1)-Bc)/(A^2) - 2/9*(2*(beta*Dx + 1)-Bc)*(11/4*Dx*beta+2-Bc)/(A^2*Dx*beta);

r=linspace(0,1,100);
x = -((B-Bc)/fi)^(1/2)*cos(nc*pi*r/L)+2*L+2/9*(((B-Bc)/fi)*(2*(beta*Dx + 1)-Bc)/(A*Dx*beta))*cos(2*nc*pi*r/L);
y = ((B-Bc)/fi)^(1/2)*(-0.6)*((Dx*beta+1-Bc)/A^2)*cos(nc*pi*r/L) + 2/9*(((B-Bc)/fi)*(2*(beta*Dx + 1)-Bc)/(A^3*Dx*beta))*((7/4)*Dx*beta + 1 - Bc)*cos(2*nc*pi*r/L) - (1/2)*((B-Bc)/fi)*(2*(beta*Dx + 1)-Bc)/A^3+2*L;

figure
plot(r,x)
title('Analytical steady state dissipative structure for X vs r')
xlabel('Distance r')
ylabel('X')

figure
plot(r,y)
title('Analytical steady state dissipative structure for Y vs r')
xlabel('Distance r')
ylabel('Y')