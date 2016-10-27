%% test_funcs
% Author: Wisam Reid 
%
%% CLEAN AND CLEAR

clc;
clear;
close all;

%% Utilities

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sets Matlab's current
% directory
%
% Adds all of the parent 
% folder's subdirectories
%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp = matlab.desktop.editor.getActive;
file = tmp.Filename;
clear tmp
cd(fileparts(file));
file(find(file=='/',1,'last'):end) = [];
addpath(genpath(file))
clear file
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Using VanDerPol.m
% 
% Numerical integration of Duffing-VanDerPol Oscillator.
% Ordinary differential equation solver using built-in
% Matlab function (ode45).
% 
% Set initial conditions for the three differential equations 
xo=[1;0;0];
% Set integration time
timePoints=5000:0.3:7500;
% Solve equations 
options=odeset('RelTol',1e-10,'AbsTol',1e-12);
[t,x]=ode45('VanDerPol',timePoints,xo);

%% Three-Dimensional Phase Space Plots
%
% Transform x3 to sin(x3) 
xnew = x; 
xnew(:,3) = sin(x(:,3));

% Generate3Dplot
plot3(xnew(:,1),xnew(:,2),xnew(:,3)); 
xlabel('x1');
ylabel('x2');
zlabel('sin(x3)')
title('Duffing-VanderPol oscillator.'); 
rotate3d on
view([-5,58]);

%% Two-DimensionalSurface of Section
%
% Two dimension Poincare surface of section on plane x3=0
% (data assumed to be stored in matrix 'xnew' generated from PSp1ot.m). 
%
% Startwithemptysurface
% clear Surface
OldDistance=0;
k=1;
for i=1:size(xnew)*[1;0]
    NewDistance = xnew(i,3);
    if (NewDistance>=0* OldDistance<0)
        % Add new point to the surface 
        TotalDistance=NewDistance-OldDistance; 
        Surface(k,:)=xnew(i-1,:)-(OldDistance/TotalDistance)*...
            (xnew(i,:) - xnew(i-1,:)); 
        k = k+1;
    end
    OldDistance=NewDistance; 
end
% Generate2Dplot 
plot(Surface(:,1),Surface(:,2),'*'); 
xlabel('x1');
ylabel('x2');
title('Poincare Surface of Section.');
