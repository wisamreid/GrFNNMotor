function [regime, fixed_pts, r_star, psi_star] = getFP(f_osc, f_input, alpha, beta1, beta2, epsilon, F)
% 
% getFP: Calculate fixed points for a given oscillator
%        under periodic forcing
% Author: Wisam Reid
% Email: wisam@ccrma.stanford.edu
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
%          regime: integer (1-5)
%       fixed_pts: Array
%          r_star: A float, steady state gain
%        psi_star: A float, steady state phase

% get the steady state amplitude and phase values
[r_star, psi_star] = getSS(f_osc, f_input, alpha, beta1, beta2, epsilon, F);

% Calculate the Jacobian
% See Kim & Large 2015 pg 4
% 
% partial derivative r_dot wrt r
dpdr = alpha + 3*beta1*r_star^2 + ...
        ((epsilon*beta2*r_star^4*(5 - 3*epsilon*r_star^2))/...
        (1 - epsilon*r_star^2)^2);
% partial derivative r_dot wrt psi
dpdP = -F*sin(psi_star);
% partial derivative psi_dot wrt r
dqdr = (F/r_star^2)*sin(psi_star);
% partial derivative psi_dot wrt psi
dqdP = -1*(F/r_star)*cos(psi_star);

Jacobian = [dpdr, dpdP; ...
            dqdr, dqdP]


end

