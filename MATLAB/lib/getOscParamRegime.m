function [regime] = getOscParamRegime(alpha, beta1, beta2, epsilon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% getOscParamRegime: 
%         Calculate which bifurcation parameter regime we 
%         are operating in
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

%% Determine parameter regime

regime = 0;


end

