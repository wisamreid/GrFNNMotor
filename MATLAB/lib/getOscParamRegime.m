function [regime] = getOscParamRegime(alpha, beta1, beta2, epsilon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% getOscParamRegime: 
%         Calculate which bifurcation parameter regime of 
%         the autonomous oscillator behavior 
% 
% Author: Wisam Reid
%  Email: wisam@ccrma.stanford.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% Function Definition:
% 
%   findOscParamRegime(alpha, beta1, beta2, epsilon)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% 
% Arguments: 
% 
%           alpha:  (float) Hopf Oscillator Parameter
%           beta1:  (float) Hopf Oscillator Parameter
%           beta2:  (float) Hopf Oscillator Parameter
%         epsilon:  (float) Hopf Oscillator Parameter
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% 
% Returns:  [regime]
% 
%          regime: integer (0-4)
%                   1: critical Hopf
%                   2: super critical Hopf
%                   3: super critical double limit cycle
%                   4: subcritical double limit cycle 
%                   0: could not be identified
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% Handle Arguments

switch nargin
    case nargin ~= 4       
        disp('Error running getOscParamRegime')
        disp('Please provide the right number of arguments')
        disp('Function Call: ')
        disp('getOscParamRegime(alpha, beta1, beta2, epsilon)')
        fprintf('\n')
        disp('type "help getOscParamRegime" for more details')
        fprintf('\n')
        return
end 

%% Define r_dot

r_dot = @(r) alpha*r + beta1*r^3 + ((epsilon*beta2*r^5)/(1 - epsilon*r^2));

%% Determine parameter regime

syms r 
% Set the derivative of r_dot wrt r equal to zero
d_r_dot_dr = alpha + ...
            (3*beta1*r^2) + ...
            ((epsilon*beta2*r^4*(5 - 3*epsilon*r^2))/((-1 + epsilon*r^2)^2)) ...
            == 0;

sol = vpasolve(d_r_dot_dr, r);

% plug the largest
local_max = r_dot(max(sol));

%% Evaluate Regime

if alpha == 0 && beta1 < 0 && beta2 == 0
    if local_max == 0    
        regime = 1;
    else
        display('Warning: critical Hopf Regime does not have')
        display('a fixed point at 0')
        regime = 0; 
    end
elseif alpha > 0 && beta1 < 0 && beta2 == 0
    if local_max > 0    
        regime = 2;
    else
        display('Warning: supercritical Hopf regime does not have')
        display('a fixed point greater than 0')
        regime = 0; 
    end
elseif alpha < 0 && beta1 > 0 && beta2 < 0
    if local_max > 0    
        regime = 3;
    elseif local_max < 0
        regime = 4;
    else
        display('Warning: Parameter regime could not be classified')
        regime = 0; 
    end
else
   display('Warning: Parameter regime could not be classified')
   regime = 0; 
end


end

