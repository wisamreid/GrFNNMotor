function poincare

%%% Compute Poincare section or First return map from a time series. 
%%% Author: Didier Gonze
%%% Created 8/3/2011
%%% Updated 22/3/2011

clear;
clc;


% ------------------------------------------------------------
% Data & parameters
% ------------------------------------------------------------

D=load('myfile.dat');   % Read data

type='max';     % 'cut', 'max', 'min', 'period', 'mean'
var=2;          % variable used for the cut/first return map (column of D)

cut=2;          % to be used with type = 'cut' (value for the cut)
dir=1;          % to be used with type = 'cut' (1=upwards;-1=downwards;0=both)

varo1=3;        % output variable 1 used with 'cut' or 'mean' (column of D)
varo2=4;        % output variable 2 used with 'cut' or 'mean' (column of D)


% ------------------------------------------------------------
% Compute First return map for max
% Return max_n as a function of max_(n-1)
% ------------------------------------------------------------

if strcmp(type,'max')

x=D(:,var);     % extract time series for variable var

xt=x(2:end-1);  % variable
xp=x(1:end-2);  % variable shifted to left
xm=x(3:end);    % variable shifted to right

kmax=find((xt>xp) & (xt>=xm));    % index of maxima

xmax=xt(kmax);       % maxima

R=[];
for i=2:length(kmax)
    R=[R; xmax(i-1) xmax(i)];
end

figure(1)
clf;
plot(R(:,1),R(:,2),'b.')
xlabel('Max_n','fontsize',18)
ylabel('Max_n_+_1','fontsize',18)
title('Return map of the maxima','fontsize',18) 

end


% ------------------------------------------------------------
% Compute First return map for min
% Return min_n as a function of min_(n-1)
% ------------------------------------------------------------

if strcmp(type,'min')

x=D(:,var);     % extract time series for variable var

xt=x(2:end-1);  % variable
xp=x(1:end-2);  % variable shifted to left
xm=x(3:end);    % variable shifted to right

kmin=find((xt<xp) & (xt<=xm));    % index of minima

xmin=xt(kmin);       % minima

R=[];
for i=2:length(kmin)
    R=[R; xmin(i-1) xmin(i)];
end

figure(1)
clf;
plot(R(:,1),R(:,2),'b.')
xlabel('Min_n','fontsize',18)
ylabel('Min_n_+_1','fontsize',18)
title('Return map of the minima','fontsize',18) 

end


% ------------------------------------------------------------
% Compute First return map for period
% Return period_n as a function of period_(n-1)
% ------------------------------------------------------------

if strcmp(type,'period')

x=D(:,var);     % extract time series for variable var
t=D(:,1);     % extract time series for time (1st column)

xt=x(2:end-1);  % variable
xp=x(1:end-2);  % variable shifted to left
xm=x(3:end);    % variable shifted to right

k=find((xt>xp) & (xt>=xm));    % index of maxima

tp=t(k);   % time of max x

R=[];
pold=tp(2)-tp(1);
for i=3:length(k)
    p=tp(i)-tp(i-1);
    R=[R; pold p];
    pold=p;
end

figure(1)
clf;
plot(R(:,1),R(:,2),'b.')
xlabel('Period_n','fontsize',18)
ylabel('Period_n_+_1','fontsize',18)
title('Return map of the period','fontsize',18) 

end


% ------------------------------------------------------------
% Compute Poicare section for cut
% Return values of all variables at cut(var) = cut
% ------------------------------------------------------------

if strcmp(type,'cut')

x=D(:,var);     % extract time series for variable var

xt=x(2:end-1);  % variable
xp=x(1:end-2);  % variable shifted to left
xm=x(3:end);    % variable shifted to right

if (dir==1) 
k=find((xp<cut) & (xt>=cut));    % index of cut points (cross upwards)
elseif (dir==-1) 
k=find((xp>cut) & (xm<=cut));    % index of cut points (cross downwards)
elseif (dir==0)
k=find(((xp<cut) & (xt>=cut)) | ((xp>cut) & (xm<=cut)));    % index of cut points (cross up- or downwards)
end

R=[];
for i=1:length(k)
    R=[R; D(k(i),:)];
end

figure(1)
clf;
plot(R(:,varo1),R(:,varo2),'b.')         % specify the two variables to be plotted
xlabel('Var (1)','fontsize',18)
ylabel('Var (2)','fontsize',18)
title(sprintf('Poincare section (cut=%g)',cut),'fontsize',18) 

end


% ------------------------------------------------------------
% Compute Poicare section for mean
% Return values of all variables at cut(var) = mean(var)
% ------------------------------------------------------------

if strcmp(type,'mean')
    
x=D(:,var);     % extract time series for variable var
cut=mean(x);    % cut at mean(x)

xt=x(2:end-1);  % variable
xp=x(1:end-2);  % variable shifted to left
xm=x(3:end);    % variable shifted to right

if (dir==1) 
k=find((xp<cut) & (xt>=cut));    % index of cut points (cross upwards)
elseif (dir==-1) 
k=find((xp>cut) & (xm<=cut));    % index of cut points (cross downwards)
elseif (dir==0)
k=find(((xp<cut) & (xt>=cut)) | ((xp>cut) & (xm<=cut)));    % index of cut points (cross up- or downwards)
end


R=[];
for i=1:length(k)
    R=[R; D(k(i),:)];
end

figure(1)
clf;
plot(R(:,varo1),R(:,varo2),'b.')         % specify the two variables to be plotted
xlabel('Var (1)','fontsize',18)
ylabel('Var (2)','fontsize',18)
title(sprintf('Poincare section (cut at mean=%g)',cut),'fontsize',18) 

end


