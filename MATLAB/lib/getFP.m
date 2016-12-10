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
%                   1: stable node
%                   2: stable spiral
%                   3: unstable node
%                   4: unstable spiral
%                   5: a saddle point
%                   0: could not be identified
%       fixed_pts: Array
%          r_star: A float, steady state gain
%        psi_star: A float, steady state phase


if nargin ~= 7
    disp('Error running getFP')
    disp('Please provide the right number of arguments')
    disp('Function Call: ')
    disp('getFP(f_osc, f_input, alpha, beta1, beta2, epsilon, F)')
    fprintf('\n')
    disp('type "help getFP" for more details')
    fprintf('\n')
    return
end 

% Displays 
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
   

% get the steady state amplitude and phase values
[r_star, psi_star] = getSS(f_osc, f_input, alpha, beta1, beta2, epsilon, F, 0);

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

% the Jacobian
J = [dpdr, dpdP; ...
            dqdr, dqdP];
        
T = trace(J); % calculate the trace
Delta = det(J); % calculate the determinant
disc = T^2 - 4*Delta; % calculate the discriminant

if Delta > 0 && disc > 0 && T < 0;
    regime = 1;
    disp(['The fixed point is a stable node'])
    fprintf('\n') 

elseif Delta > 0 && disc < 0 && T < 0;
    regime = 2;
    disp(['The fixed point is a stable spiral'])
    fprintf('\n')

elseif Delta > 0 && disc > 0 && T > 0;
    regime = 3;
    disp(['The fixed point is a unstable node'])
    fprintf('\n') 

elseif Delta > 0 && disc < 0 && T > 0;
    regime = 4;
    disp(['The fixed point is a unstable spiral'])
    fprintf('\n') 

elseif Delta < 0;
    regime = 5;
    disp(['The fixed point is a saddle point'])
    disp(['See Strogatz 1994'])
    fprintf('\n')

else Delta < 0;
    regime = 0;
    disp(['The fixed point could not be identified'])
    fprintf('\n')
end 

%  add fixed points later
fixed_pts = ['fixed points have not been implimented yet'];

end

