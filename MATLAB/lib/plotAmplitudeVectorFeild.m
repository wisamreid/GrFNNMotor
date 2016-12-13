function [ ] = plotAmplitudeVectorFeild(alpha, beta1, beta2, epsilon)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% plotAmplitudeVectorFeild: 
%         Plots the amplitude vector feild for a given
%         oscillator
% 
% Author: Wisam Reid
%  Email: wisam@ccrma.stanford.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% Function Definition:
% 
%   plotAmplitudeVectorFeild(alpha, beta1, beta2, epsilon)
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
% Returns:  N/A
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Handle Arguments

switch nargin
    case nargin ~= 4       
        disp('Error running plotAmplitudeVectorFeild')
        disp('Please provide the right number of arguments')
        disp('Function Call: ')
        disp('plotAmplitudeVectorFeild(alpha, beta1, beta2, epsilon)')
        fprintf('\n')
        disp('type "help plotAmplitudeVectorFeild" for more details')
        fprintf('\n')
        return
end 

%% Define r_dot

r_dot = @(r) alpha.*r + beta1.*r.^3 + ((epsilon*beta2.*r.^5)/(1 - epsilon.*r.^2));

r_dot(0)

% r values
r = 0:0.01:1;

%% Plot

regime = getOscParamRegime(alpha, beta1, beta2, epsilon); % Check for a known regime

amp = r_dot(r);

figure(1)
if regime == 1
    plot(r,amp)
    hold on 
    plot([0 max(r)],[0 0],'k-')
    ylim([min(amp),0.2])
elseif regime == 2
    plot(r,amp)
    hold on 
    plot([0 max(r)],[0 0],'k-')
    ylim([min(amp),max(amp)+ 0.1])
elseif regime == 3
    plot(r,amp)
    hold on 
    plot([0 max(r)],[0 0],'k-')
    ylim([min(amp),max(amp)+ 0.1])
else
    display('Warning: The oscillator parameter regime could not be classified')
end

end

