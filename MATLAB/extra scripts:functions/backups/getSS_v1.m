function [r_star, psi_star] = getSS(f_osc, f_input, alpha, beta1, beta2, epsilon, F, display_flag)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
%  getSS: Calculate steady state amplitude and phase 
%         for a given oscillator under periodic forcing
% 
% Author: Wisam Reid
%  Email: wisam@ccrma.stanford.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% Function Definition:
% 
%   getSS(f_osc, f_input, alpha, beta1, beta2, epsilon, F, display_flag)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% 
% Arguments: 
% 
%           f_osc:  (float) Oscillator Frequency
%         f_input:  (float) Forcing Input Frequency
%           alpha:  (float) Hopf Oscillator Parameter
%           beta1:  (float) Hopf Oscillator Parameter
%           beta2:  (float) Hopf Oscillator Parameter
%         epsilon:  (float) Hopf Oscillator Parameter
%               F:  (float) Forcing Amplitude
%    display_flag:  [optional] (boolean) print to console
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% 
% Returns:  [r_star, psi_star,regime]
% 
%          r_star: A float, steady state gain
%        psi_star: A float, steady state phase
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% Handle Arguments

switch nargin
    case nargin == 7
        display_flag = 1;
    case nargin ~= 8 ||  nargin ~= 7       
        disp('Error running getSS')
        disp('Please provide the right number of arguments')
        disp('Function Call: ')
        disp('getSS(f_osc, f_input, alpha, beta1, beta2, epsilon, F, print_flag)')
        fprintf('\n')
        disp('type "help getSS" for more details')
        fprintf('\n')
        return
end 

% Check input arguments
if F <= 0
  error('F > 0')
end

%% Displays

if display_flag
    % Display
    disp(['Running getFP:'])
    fprintf('\n')
    disp(['The Oscillator Frequency is: ',num2str(f_osc), ' Hz'])
    disp(['The Input Frequency is: ',num2str(f_input), ' Hz'])
    fprintf('\n')
    disp(['alpha is: ',num2str(alpha)])
    disp(['beta1 is: ',num2str(beta1)])
    disp(['beta2 is: ',num2str(beta2)])
    disp(['epsilon is: ',num2str(epsilon)])
    disp(['F is: ',num2str(F)])
    fprintf('\n')
end

%% Determine parameter regime

% % First, we need to know which bifurcation parameter regime we 
% % are operating in
% regime = getOscParamRegime(alpha, beta1, beta2, epsilon);

%%
% Next, we can find the boundary between stable nodes and stable spirals by solving:
%     
%     discriminant = 0
%     r_dot = 0 
%     psi_dot = 0 
%     
% simultaneously where T and   are the trace and determinant of the 
% Jacobian matrix evaluated at a fixed point (see Section 2). 
% We find that the boundary is at

%%
% calculate Omega
Omega = 2*pi*(f_osc - f_input);

syms r Psi
% z_dot = [alpha*r + beta1*r^3 + ((epsilon*beta2*r^5)/(1 - epsilon*r^2)) + F*cos(Psi) == 0, ...
%         Omega - (F/r)*sin(Psi) == 0];

z_sol = vpasolve([alpha*r + beta1*r^3 + ...
             ((epsilon*beta2*r^5)/(1 - epsilon*r^2)) + F*cos(Psi) == 0, ...
             Omega - (F/r)*sin(Psi) == 0], [r, Psi],'random',true);

% % % function F = root2d(x);
% % %     
% % % F(1) = alpha*r + beta1*r^3 + ...            
% % %        ((epsilon*beta2*r^5)/(1 - epsilon*r^2)) + F*cos(Psi);
% % % F(2) = Omega - (F/r)*sin(Psi);
% % % 
% % % tspan = 0:0.1:10;  %you might want to adapt the range for t
% % % [T,Z] = ode45(@myfunc, tspan, [x0 y0]); %myfun is your DE function
% % % x = Z(1,:);
% % % y = Z(2,:);

% r_star = double(vpa(z_sol.r))
% psi_star = double(eval(mod(vpa(z_sol.Psi), 2*pi))) % wrap the phase

r_star = double(z_sol.r)
psi_star = double(eval(mod(z_sol.Psi, 2*pi))) % wrap the phase

%% Now we handle each quadrant correctly

% first quadrant
% do nothing

% second quadrant
if psi_star > pi/2 && psi_star < pi
    % do nothing
end

% third quadrant 
if psi_star > pi && psi_star < 3*pi/2
    psi_star = psi_star - pi;
    r_star = -1*r_star;
end

% fourth quadrant 
if psi_star > 3*pi/2
    psi_star = mod(psi_star + pi, 2*pi);
    r_star = -1*r_star;
end

% edge cases
% if psi_star == 0 || psi_star == pi
%    r_star = -1*r_star; 
% end

% % print
% r_star
% psi_star

end

