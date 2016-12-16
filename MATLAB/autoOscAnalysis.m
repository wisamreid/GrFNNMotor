%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% Autonomous Oscillator Parameter Analysis: 
%         Plots the amplitude vector feild for ranging 
%         oscillator parameters
% 
% Author: Wisam Reid
%  Email: wisam@ccrma.stanford.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%% CLEAN AND CLEAR

close all 
clear
clc

%% Plots critical hopf sweep alpha

plot_number = 1;

alpha = 0;
beta1 = 0;
beta2 = 0;
epsilon = 0; 
F = 0;
f_osc = 1;
f_input = 0;

for beta1 = 1:10
    fig1 = plotAmplitudeVectorFeild(alpha, -1*beta1, beta2, epsilon,...
                                    F, f_osc, f_input, plot_number);
end
plot([0 1],[0 0],'k--','linewidth',2)
hold on
plotFPs(0, -1*beta1, 0, 0);
axis([0 1 -0.01 0.002])
title('Autonomous Critical Oscillator: $$\beta_1$$ sweep: $$\alpha = 0$$, $$\beta_2 = 0$$, $$\epsilon = 0$$',...
      'Interpreter', 'latex')
legend({'$$\beta_1 = 1$$','$$\beta_1 = 2$$','$$\beta_1 = 3$$', ... 
        '$$\beta_1 = 4$$','$$\beta_1 = 5$$','$$\beta_1 = 6$$', ...
        '$$\beta_1 = 7$$','$$\beta_1 = 8$$','$$\beta_1 = 9$$', '$$\beta_1 = 10$$'},...
        'Location','southwest', 'Interpreter', 'latex');
grid on    
hold off

%% Plots supercritical hopf sweep alpha

plot_number = 2;

alpha = 0;
beta1 = 0;
beta2 = 0;
epsilon = 0; 
F = 0;
f_osc = 1;
f_input = 0;


for beta1 = 1:10
    fig1 = plotAmplitudeVectorFeild(alpha, -1*beta1, beta2, epsilon,...
                                    F, f_osc, f_input, plot_number);
end
plot([0 1],[0 0],'k--','linewidth',2)
hold on
plotFPs(0, -1*beta1, 0, 0);
axis([0 1 -0.01 0.002])
title('Autonomous Critical Oscillator: $$\beta_1$$ sweep: $$\alpha = 0$$, $$\beta_2 = 0$$, $$\epsilon = 0$$',...
      'Interpreter', 'latex')
legend({'$$\beta_1 = 1$$','$$\beta_1 = 2$$','$$\beta_1 = 3$$', ... 
        '$$\beta_1 = 4$$','$$\beta_1 = 5$$','$$\beta_1 = 6$$', ...
        '$$\beta_1 = 7$$','$$\beta_1 = 8$$','$$\beta_1 = 9$$', '$$\beta_1 = 10$$'},...
        'Location','southwest', 'Interpreter', 'latex');
grid on    
hold off
