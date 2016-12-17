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

% oscillator params
alpha = 0;
beta1 = 0;
beta2 = 0;
epsilon = 0; 
F = 0;
f_osc = 1;
f_input = 0;

sweep = 1:10;
numSteps = length(sweep);

C = distinguishable_colors(numSteps);

i = 1;
for beta1 = sweep
    
    fig1 = plotAmplitudeVectorFeild(alpha, -1*beta1, beta2, epsilon,...
                                    F, f_osc, f_input, plot_number, C(i,:));
    i = i + 1 ;
    
end
plot([0 1],[0 0],'k-','linewidth',2)
hold on
plotFPs(alpha, -1*beta1, beta2, epsilon,...
        F, f_osc, f_input, plot_number);
axis([0 0.5 -0.001 0.00010])
title({'Amplitude Vector Field of an Autonomous Critical-Hopf Oscillator: $$\beta_1$$ sweep';...
        '$$\qquad \qquad \qquad \qquad \qquad \qquad \alpha = 0$$, $$\beta_2 = 0$$, $$\epsilon = 0$$'}, ...
        'Interpreter', 'latex')
legend({'$$\beta_1 = -1$$','$$\beta_1 = -2$$','$$\beta_1 = -3$$', ... 
        '$$\beta_1 = -4$$','$$\beta_1 = -5$$','$$\beta_1 = -6$$', ...
        '$$\beta_1 = -7$$','$$\beta_1 = -8$$','$$\beta_1 = -9$$', '$$\beta_1 = -10$$'},...
        'Location','southwest', 'Interpreter', 'latex');
xlabel('$$r$$', 'Interpreter', 'latex')
ylabel('$$\dot{r}$$', 'Interpreter', 'latex')
grid on    
hold off

%% Plots supercritical hopf sweep alpha hold beta1

plot_number = 2;

% oscillator params
alpha = 1;
beta1 = -10;
beta2 = 0;
epsilon = 0; 
F = 0;
f_osc = 1;
f_input = 0;

sweep = 1:10;
numSteps = length(sweep);

C = distinguishable_colors(numSteps);

i = 1;
for alpha = sweep
    
    fig2 = plotAmplitudeVectorFeild(alpha, beta1, beta2, epsilon,...
                                    F, f_osc, f_input, plot_number, C(i,:));
    i = i + 1;
    
end
plot([0 1],[0 0],'k-','linewidth',2)
hold on
for alpha = sweep
    plotFPs(alpha, beta1, beta2, epsilon,...
            F, f_osc, f_input, plot_number);
end
axis([0 1 -0.0015 0.004])
% title('Autonomous Supercritical Oscillator: $$\alpha$$ sweep: $$\beta_1 = -10$$, $$\beta_2 = 0$$, $$\epsilon = 0$$',...
%       'Interpreter', 'latex')
title({'Amplitude Vector Field of an Autonomous Supercritical Oscillator: $$\alpha$$ sweep';...
        '$$\qquad \qquad \qquad \qquad \qquad \qquad \beta_1 = -10$$, $$\beta_2 = 0$$, $$\epsilon = 0$$'}, ...
        'Interpreter', 'latex')
legend({'$$\alpha = 1$$','$$\alpha = 2$$','$$\alpha = 3$$', ... 
        '$$\alpha = 4$$','$$\alpha = 5$$','$$\alpha = 6$$', ...
        '$$\alpha = 7$$','$$\alpha = 8$$','$$\alpha = 9$$', '$$\alpha = 10$$'},...
        'Location','northwest', 'Interpreter', 'latex');
xlabel('$$r$$', 'Interpreter', 'latex')
ylabel('$$\dot{r}$$', 'Interpreter', 'latex')
grid on    
hold off

%% Plots supercritical hopf sweep beta1 hold alpha

plot_number = 3;

% oscillator params
alpha = 10;
beta1 = 0;
beta2 = 0;
epsilon = 0; 
F = 0;
f_osc = 1;
f_input = 0;

sweep = 1:10;
numSteps = length(sweep);

C = distinguishable_colors(numSteps);

i = 1;
for beta1 = sweep
    
    fig3 = plotAmplitudeVectorFeild(alpha, -1*beta1, beta2, epsilon,...
                                    F, f_osc, f_input, plot_number, C(i,:));
    i = i + 1;
    
end
plot([0 3.5],[0 0],'k-','linewidth',2)
hold on
for beta1 = 1:10
    plotFPs(alpha,-1*beta1, beta2, epsilon,...
            F, f_osc, f_input, plot_number);
end
axis([0 3.5 -0.004 0.014])
% axis([0 1 -0.0001 0.0005])
title({'Amplitude Vector Field of an Autonomous Supercritical Oscillator: $$\beta_1$$ sweep';...
        '$$\qquad \qquad \qquad \qquad \qquad \qquad \alpha = 10$$, $$\beta_2 = 0$$, $$\epsilon = 0$$'}, ...
        'Interpreter', 'latex')
legend({'$$\beta_1 = -1$$','$$\beta_1 = -2$$','$$\beta_1 = -3$$', ... 
        '$$\beta_1 = -4$$','$$\beta_1 = -5$$','$$\beta_1 = -6$$', ...
        '$$\beta_1 = -7$$','$$\beta_1 = -8$$','$$\beta_1 = -9$$', '$$\beta_1 = -10$$'},...
        'Location','northeast', 'Interpreter', 'latex');
xlabel('$$r$$', 'Interpreter', 'latex')
ylabel('$$\dot{r}$$', 'Interpreter', 'latex')
grid on    
hold off

%% Plots supercritical hopf sweep beta1 and alpha
% Not for the current research report 

% plot_number = 4;
% 
% % oscillator params
% alpha = 10;
% beta1 = 0;
% beta2 = 0;
% epsilon = 0; 
% F = 0;
% f_osc = 1;
% f_input = 0;
% 
% sweep = 1:10;
% numSteps = length(sweep);
% 
% C = distinguishable_colors(numSteps);
% 
% i = 1;
% for beta1 = sweep
%     alpha = 11 - beta1;
%     fig4 = plotAmplitudeVectorFeild(alpha, -1*beta1, beta2, epsilon,...
%                                     F, f_osc, f_input, plot_number, C(i,:));
%     i = i + 1;
% end
% plot([0 3.5],[0 0],'k-','linewidth',2)
% hold on
% for beta1 = sweep
%     alpha = 11 - beta1;
%     plotFPs(alpha,-1*beta1, beta2, epsilon,...
%             F, f_osc, f_input, plot_number);
% end
% axis([0 3.5 -0.0015 0.014])
% title({'Amplitude Vector Field of an Autonomous Supercritical Oscillator: $$\alpha$$ and $$\beta_1$$ sweep';...
%         '$$\qquad \qquad \qquad \qquad \qquad \qquad \qquad \qquad \beta_2 = 0$$, $$\epsilon = 0$$'}, ...
%         'Interpreter', 'latex')
% legend({'$$\alpha = 1,~ \beta_1 = 10$$','$$\alpha = 2,~ \beta_1 = 9$$',...
%         '$$\alpha = 3,~ \beta_1 = 8$$', '$$\alpha = 4,~ \beta_1 = 7$$',...
%         '$$\alpha = 5,~ \beta_1 = 6$$','$$\alpha = 6,~ \beta_1 = 5$$', ...
%         '$$\alpha = 7,~ \beta_1 = 4$$','$$\alpha = 8,~ \beta_1 = 3$$', ...
%         '$$\alpha = 9,~ \beta_1 = 2$$', '$$\alpha = 10,~ \beta_1 = 1$$'},...
%         'Location','northeast', 'Interpreter', 'latex');
% grid on    
% hold off

%% Plots supercritical DLC sweep alpha

plot_number = 5;

% oscillator params
alpha = 0;
beta1 = 10;
beta2 = -1;
epsilon = 1; 
F = 0;
f_osc = 1;
f_input = 0;


sweep = 0.2:0.2:2;
numSteps = length(sweep);

C = distinguishable_colors(numSteps);

i = 1;
for alpha = sweep
    
    fig5 = plotAmplitudeVectorFeild(-1*alpha, beta1, beta2, epsilon,...
                                    F, f_osc, f_input, plot_number, C(i,:));
    i = i + 1;
end
plot([0 1],[0 0],'k-','linewidth',2)
hold on
for alpha = sweep
    
    plotFPs(-1*alpha, beta1, beta2, epsilon,...
            F, f_osc, f_input, plot_number);
end
axis([0 1 -0.005 0.005])
title({'Amplitude Vector Field of an Autonomous Supercritical DLC Oscillator: $$\alpha$$ sweep';...
        '$$\qquad \qquad \qquad \qquad \qquad \qquad \beta_1 = 10$$, $$\beta_2 = -1$$, $$\epsilon = 1$$'}, ...
        'Interpreter', 'latex')
legend({'$$\alpha = -0.2$$','$$\alpha = -0.4$$','$$\alpha = -0.6$$', '$$\alpha = -0.8$$',...
        '$$\alpha = -1$$','$$\alpha = -1.2$$','$$\alpha = -1.4$$','$$\alpha = -1.6$$', ...
        '$$\alpha = -1.8$$', '$$\alpha = -2$$'},...
        'Location','southwest', 'Interpreter', 'latex');
xlabel('$$r$$', 'Interpreter', 'latex')
ylabel('$$\dot{r}$$', 'Interpreter', 'latex')
grid on    
hold off

%% Plots supercritical DLC sweep beta1 

plot_number = 6;

% oscillator params
alpha = -2;
beta1 = 10;
beta2 = -1;
epsilon = 1; 
F = 0;
f_osc = 1;
f_input = 0;

sweep = 9.2:0.2:11;
numSteps = length(sweep);

C = distinguishable_colors(numSteps);

i = 1;
for beta1 = sweep
    
    fig6 = plotAmplitudeVectorFeild(alpha, beta1, beta2, epsilon,...
                                    F, f_osc, f_input, plot_number, C(i,:));                           
    i = i + 1;
    
end
plot([0 1],[0 0],'k-','linewidth',2)
hold on
for beta1 = sweep
    
    plotFPs(alpha, beta1, beta2, epsilon,...
            F, f_osc, f_input, plot_number);
end
axis([0 1 -0.005 0.005])
title({'Amplitude Vector Field of an Autonomous Supercritical DLC Oscillator: $$\beta_1$$ sweep';...
        '$$\qquad \qquad \qquad \qquad \qquad \qquad \alpha = -2$$, $$\beta_2 = -1$$, $$\epsilon = 1$$'}, ...
        'Interpreter', 'latex')
legend({'$$\beta_1 = 9.2$$','$$\beta_1 = 9.4$$','$$\beta_1 = 9.6$$', '$$\beta_1 = 9.8$$',...
        '$$\beta_1 = 10$$','$$\beta_1 = 10.2$$','$$\beta_1 = 10.4$$','$$\beta_1 = 10.6$$', ...
        '$$\beta_1 = 10.8$$', '$$\beta_1 = 11$$'},...
        'Location','southwest', 'Interpreter', 'latex');
xlabel('$$r$$', 'Interpreter', 'latex')
ylabel('$$\dot{r}$$', 'Interpreter', 'latex')
grid on    
hold off

%% Plots supercritical DLC sweep beta2 

plot_number = 7;

% oscillator params
alpha = -2;
beta1 = 10;
beta2 = -1;
epsilon = 1; 
F = 0;
f_osc = 1;
f_input = 0;

sweep = 0.5:0.1:1.4;
numSteps = length(sweep);

C = distinguishable_colors(numSteps);

i = 1;
for beta2 = sweep
    
    fig5 = plotAmplitudeVectorFeild(alpha, beta1, -1*beta2, epsilon,...
                                    F, f_osc, f_input, plot_number, C(i,:));
    i = i + 1;
    
end
plot([0 1],[0 0],'k-','linewidth',2)
hold on
for beta2 = sweep
    
    plotFPs(alpha, beta1, -1*beta2, epsilon,...
            F, f_osc, f_input, plot_number);
end
axis([0 1 -0.005 0.005])
title({'Amplitude Vector Field of an Autonomous Supercritical DLC Oscillator: $$\beta_2$$ sweep';...
        '$$\qquad \qquad \qquad \qquad \qquad \qquad \alpha = -2$$, $$\beta_1 = 10$$, $$\epsilon = 1$$'}, ...
        'Interpreter', 'latex')
legend({'$$\beta_2 = -0.5$$','$$\beta_2 = -0.6$$','$$\beta_2 = -0.7$$', '$$\beta_2 = -0.8$$',...
        '$$\beta_2 = -0.9$$','$$\beta_2 = -1$$','$$\beta_2 = -1.1$$','$$\beta_2 = -1.2$$', ...
        '$$\beta_2 = -1.3$$','$$\beta_2 = -1.4$$'},...
        'Location','southwest', 'Interpreter', 'latex');
xlabel('$$r$$', 'Interpreter', 'latex')
ylabel('$$\dot{r}$$', 'Interpreter', 'latex')
grid on    
hold off

%% Plots subcritical DLC sweep alpha

plot_number = 8;

% oscillator params
alpha = -6;
beta1 = 8;
beta2 = -1;
epsilon = 1; 
F = 0;
f_osc = 1;
f_input = 0;

sweep = 6:0.2:7.8;
numSteps = length(sweep);

C = distinguishable_colors(numSteps);

i = 1;
for alpha = sweep
    
    fig5 = plotAmplitudeVectorFeild(-1*alpha, beta1, beta2, epsilon,...
                                    F, f_osc, f_input, plot_number, C(i,:));
    i = i + 1;
end
plot([0 1],[0 0],'k-','linewidth',2)
hold on
for alpha = sweep
    
    plotFPs(-1*alpha, beta1, beta2, epsilon,...
            F, f_osc, f_input, plot_number);
end
axis([0 1 -0.005 0.001])
title({'Amplitude Vector Field of an Autonomous Subcritical DLC Oscillator: $$\alpha$$ sweep';...
        '$$\qquad \qquad \qquad \qquad \qquad \qquad \beta_1 = 8$$, $$\beta_2 = -1$$, $$\epsilon = 1$$'}, ...
        'Interpreter', 'latex')
legend({'$$\alpha = -6$$','$$\alpha = -6.2$$','$$\alpha = -6.4$$', '$$\alpha = -6.6$$',...
        '$$\alpha = -6.8$$','$$\alpha = -7.0$$','$$\alpha = -7.2$$','$$\alpha = -7.4$$', ...
        '$$\alpha = -7.6$$', '$$\alpha = -7.8$$'},...
        'Location','southwest', 'Interpreter', 'latex');
xlabel('$$r$$', 'Interpreter', 'latex')
ylabel('$$\dot{r}$$', 'Interpreter', 'latex')
grid on    
hold off

%% Plots subcritical DLC sweep beta1

plot_number = 9;

% oscillator params
alpha = -6;
beta1 = 8;
beta2 = -1;
epsilon = 1; 
F = 0;
f_osc = 1;
f_input = 0;

sweep = 8:0.2:9.8;
numSteps = length(sweep);

C = distinguishable_colors(numSteps);

i = 1;
for beta1 = sweep
    
    fig5 = plotAmplitudeVectorFeild(alpha, beta1, beta2, epsilon,...
                                    F, f_osc, f_input, plot_number, C(i,:));
    i = i + 1;
end
plot([0 1],[0 0],'k-','linewidth',2)
hold on
for beta1 = sweep
    
    plotFPs(alpha, beta1, beta2, epsilon,...
            F, f_osc, f_input, plot_number);
end
axis([0 1 -0.005 0.001])
title({'Amplitude Vector Field of an Autonomous Subcritical DLC Oscillator: $$\beta_1$$ sweep';...
        '$$\qquad \qquad \qquad \qquad \qquad \qquad \alpha = -6$$, $$\beta_2 = -1$$, $$\epsilon = 1$$'}, ...
        'Interpreter', 'latex')
legend({'$$\beta_1 = 8$$','$$\beta_1 = 8.2$$','$$\beta_1 = 8.4$$', '$$\beta_1 = 6.6$$',...
        '$$\beta_1 = 8.8$$','$$\beta_1 = 9.0$$','$$\beta_1 = 9.2$$','$$\beta_1 = 7.4$$', ...
        '$$\beta_1 = 9.6$$', '$$\beta_1 = 9.8$$'},...
        'Location','southwest', 'Interpreter', 'latex');
xlabel('$$r$$', 'Interpreter', 'latex')
ylabel('$$\dot{r}$$', 'Interpreter', 'latex')
grid on    
hold off

%% Plots subcritical DLC sweep beta2

plot_number = 10;

% oscillator params
alpha = -6;
beta1 = 8;
beta2 = -0.2;
epsilon = 1; 
F = 0;
f_osc = 1;
f_input = 0;

sweep = 0.4:0.1:1.4;
numSteps = length(sweep);

C = distinguishable_colors(numSteps);

i = 1;
for beta2 = sweep
    
    fig5 = plotAmplitudeVectorFeild(alpha, beta1, -1*beta2, epsilon,...
                                    F, f_osc, f_input, plot_number, C(i,:));
    i = i + 1;
end
plot([0 1],[0 0],'k-','linewidth',2)
hold on
for beta2 = sweep
    
    plotFPs(alpha, beta1, -1*beta2, epsilon,...
            F, f_osc, f_input, plot_number);
end
axis([0 1 -0.005 0.001])
title({'Amplitude Vector Field of an Autonomous Subcritical DLC Oscillator: $$\beta_2$$ sweep';...
        '$$\qquad \qquad \qquad \qquad \qquad \qquad \alpha = -6$$, $$\beta_1 = 8$$, $$\epsilon = 1$$'}, ...
        'Interpreter', 'latex')
legend({'$$\beta_2 = -0.4$$','$$\beta_2 = -0.5$$','$$\beta_2 = -0.6$$', '$$\beta_2 = -0.8$$',...
        '$$\beta_2 = -0.9$$','$$\beta_2 = -1$$','$$\beta_2 = -1.1$$','$$\beta_2 = -1.2$$', ...
        '$$\beta_2 = -1.3$$', '$$\beta_2 = -1.4$$'},...
        'Location','southwest', 'Interpreter', 'latex');
xlabel('$$r$$', 'Interpreter', 'latex')
ylabel('$$\dot{r}$$', 'Interpreter', 'latex')
grid on    
hold off