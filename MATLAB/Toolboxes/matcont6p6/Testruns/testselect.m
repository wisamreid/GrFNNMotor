function testselect

OPTIONS = [];
hls = adapt2;
[t,y] = ode45(hls{2},[0 300],[0.3 0.5 -0.1],OPTIONS,1,0.8);
x0 = y(end,:);
[t,y] = ode45(hls{2},[0 10],x0,OPTIONS,1,0.8);

figure
plot(y(:,1),y(:,2))

[x,v,s,h,f] = selectcycle(t,y,1e-2,20,[1;0.8]);
par=s(end).data.parametervalues,
[x0,v0]=init_LC_LC(@adapt2,x,v,s(end),par,2,20,4);
[xlcc,vlcc,slcc,hlcc,flcc]=cont(@limitcycle,x0,v0,OPTIONS);

figure
axes
plotcycle(xlcc,vlcc,slcc,[size(xlcc,1) 1 2]);