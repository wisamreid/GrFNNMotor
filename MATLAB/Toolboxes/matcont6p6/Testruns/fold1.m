p=[0.11047;0.1];ap1=[1];
[x0,v0]=init_EP_EP(@MLfast,[0.047222;0.32564],p,ap1);
opt=contset;
opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxNumPoints',65);
opt=contset(opt,'MinStepsize',0.00001);
opt=contset(opt,'MaxStepsize',0.01);
opt=contset(opt,'Backward',1);

[x,v,s,h,f]=cont(@equilibrium,x0,[],opt);

x1=x(1:2,s(2).index);p=[x(end,s(2).index);0.1];
[x0,v0]=init_H_LC(@MLfast,x1,p,ap1,0.0001,30,4);

opt=contset;
opt=contset(opt,'IgnoreSingularity',1);
opt=contset(opt,'Singularities',1);
opt=contset(opt,'MaxNumPoints',50);
[x2,v2,s2,h2,f2]=cont(@limitcycle,x0,v0,opt);

plotcycle(x2,v2,s2,[1 2]);
