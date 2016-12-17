%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% THIS SCRIPT IS INCOMPLETE:
% 
%       Only the critical hopf sweeping F
%       case is currently implimented
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%  
% Forced Oscillator Parameter Analysis: 
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

%% Plots critical hopf sweep Forcing

plot_number = 1;

% oscillator params
alpha = 0;
beta1 = -1;
beta2 = 0;
epsilon = 0; 
F = 0.1;
f_osc = 1;
f_input = 1;

sweep = 0.1:0.1:2;
numSteps = length(sweep);

C = distinguishable_colors(numSteps);

i = 1;
for F = sweep
    
    fig1 = plotAmplitudeVectorField(alpha, beta1, beta2, epsilon,...
                                    F, f_osc, f_input, plot_number, C(i,:));
    i = i + 1 ;
    
end
plot([0 2],[0 0],'k-','linewidth',2)
hold on
for F = sweep
    plotFPs(alpha, beta1, beta2, epsilon,...
            F, f_osc, f_input, plot_number);
end
axis([0 1.3 -0.003 0.0015])
title({'Amplitude Vector Field of an Forced Critical-Hopf Oscillator: $$F$$ sweep';...
        '$$\qquad \qquad \qquad \Omega = 0$$, $$\alpha = 0$$, $$\beta_1 = -100$$, $$\beta_2 = 0$$, $$\epsilon = 0$$'}, ...
        'Interpreter', 'latex')
legend({'$$F = 0.1$$','$$F = 0.2$$','$$F = 0.3$$', ... 
        '$$F = 0.4$$','$$F = 0.5$$','$$F = 0.6$$', ...
        '$$F = 0.7$$','$$F = 0.8$$','$$F = 0.9$$', '$$F = 1$$'},...
        'Location','southwest', 'Interpreter', 'latex');
xlabel('$$r$$', 'Interpreter', 'latex')
ylabel('$$\dot{r}$$', 'Interpreter', 'latex')
grid on    
hold off

%% Plots supercritical hopf sweep Forcing
% 
% plot_number = 2;
% 
% % oscillator params
% alpha = 1;
% beta1 = -100;
% beta2 = 0;
% epsilon = 0; 
% F = 0.1;
% f_osc = 1;
% f_input = 1;
% 
% % [regime, r_local_max, drdt_local_max] = getOscParamRegime(alpha, beta1, beta2, epsilon)
% 
% 
% sweep = 0.02:0.01:0.2;
% numSteps = length(sweep);
% 
% C = distinguishable_colors(numSteps);
% 
% i = 1;
% for F = sweep
%     F
%     fig1 = plotAmplitudeVectorField(alpha, beta1, beta2, epsilon,...
%                                     F, f_osc, f_input, plot_number, C(i,:));
%     i = i + 1 ;
%     
% end
% plot([0 2],[0 0],'k-','linewidth',2)
% hold on
% for F = sweep
%     plotFPs(alpha, beta1, beta2, epsilon,...
%             F, f_osc, f_input, plot_number);
% end
% axis([0 1.3 -0.003 0.0015])
% title({'Amplitude Vector Field of an Forced Critical-Hopf Oscillator: $$F$$ sweep';...
%         '$$\qquad \qquad \qquad \Omega = 0$$, $$\alpha = 0$$, $$\beta_1 = -100$$, $$\beta_2 = 0$$, $$\epsilon = 0$$'}, ...
%         'Interpreter', 'latex')
% legend({'$$F = 0.1$$','$$F = 0.2$$','$$F = 0.3$$', ... 
%         '$$F = 0.4$$','$$F = 0.5$$','$$F = 0.6$$', ...
%         '$$F = 0.7$$','$$F = 0.8$$','$$F = 0.9$$', '$$F = 1$$'},...
%         'Location','southwest', 'Interpreter', 'latex');
% xlabel('$$r$$', 'Interpreter', 'latex')
% ylabel('$$\dot{r}$$', 'Interpreter', 'latex')
% grid on    
% hold off