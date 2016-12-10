function [ r_star, psi_star ] = getSS(f_osc, f_input, alpha, beta1, beta2, epsilon, F)
% 
% findSS: Calculate steady state amplitude and phase 
% 
% Arguments: 
% 
%           f_osc:  Oscillator Frequency
%         f_input:  Forcing Input Frequency
%           alpha: (Hopf Oscillator Parameter)
%           beta1: (Hopf Oscillator Parameter)
%           beta2: (Hopf Oscillator Parameter)
%         epsilon: (Hopf Oscillator Parameter)
%               F: Forcing Amplitude
% 
% Returns:
% 
% r_star: A float, steady state gain
% psi_star: A float, steady state phase

% calculate Omega
Omega = 2*pi*(f_osc - f_input);

syms r Psi
r_dot = alpha*r + beta1*r^3 + ((epsilon*beta2*r^5)/(1 - epsilon*r^2)) + F*cos(Psi) == 0;
psi_dot = Omega - (F/r)*sin(Psi) == 0;

% sol = solve([r_dot, psi_dot], [r, Psi]);
sol = vpasolve([r_dot, psi_dot], [r, Psi]);

r_star = sol.r;
psi_star = eval(mod(sol.Psi, 2*pi));

% first quadrant
% do nothing

% second quadrant
if psi_star > pi/2 && psi_star < pi
    % do nothing
end

% third quadrant 
if psi_star > pi && psi_star < 3*pi/2
    r_star = -1*r_star;
    psi_star = psi_star - pi;
end

% fourth quadrant 
if psi_star > 3*pi/2
    psi_star = mod(psi_star + pi, 2*pi);
    r_star = -1*r_star;
end

if psi_star == 0 || psi_star == pi
   r_star = -1*r_star; 
end

% r_star
% psi_star

end

