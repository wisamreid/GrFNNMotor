function [plot_obj, FPs] = plotFPs(alpha, beta1, beta2, epsilon, F, f_osc, f_input, figure_number)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% plotFPs: 
%         Plots the fixed points for a given
%         oscillator
%           
%   Note: * This function returns the plot as an object, allowing for 
%           manipulation of the figure outside of the function
%         * This function does not close the figure allowing for plotting
%           and parameter analysis over for loops
% 
% Colors:
%           Orange: Stable Node
%           Yellow: Stable Spiral
%            Green: Unstable Node
%             Cyan: Unstable Spiral
%          Magenta: Saddle Point
% 
% Author: Wisam Reid
%  Email: wisam@ccrma.stanford.edu
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% 
% Function Definition:
% 
%   plotAmplitudeVectorFeild(alpha, beta1, beta2, epsilon, F, f_osc, f_input, figure_number)
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% 
% Arguments: 
% 
%           alpha:  (float) Hopf Oscillator Parameter
%           beta1:  (float) Hopf Oscillator Parameter
%           beta2:  (float) Hopf Oscillator Parameter
%         epsilon:  (float) Hopf Oscillator Parameter
%               F:  [Optional] (float) 0 means autonomous 
%           f_osc:  [Optional] (float) oscillator frequency
%         f_input:  [Optional] (float) input frequency
%   figure_number:  [Optional] (int) number to be used for the plot
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% 
% Returns:  [plot_obj, FPs]
%   
%        plot_obj:  Figure (figure_number) with properties
%             FPs:  (float) Array of fixed points (r*)
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Example Calls:  
% 
% An autonomous critical Hopf oscillator (See Kim & Large 2015 Fig 1):
% 
%       fig = plotAmplitudeVectorFeild(0, -100, 0, 0)
% 
% A forced supercritical Hopf oscillator (See Kim & Large 2015 Fig 4):
%       
%       fig = plotAmplitudeVectorFeild(1, -100, 0, 0, 0.02, 1, 0.5)
% 
% A forced supercritical DLC oscillator using plot number 3:
%   (See Kim & Large 2015 Fig 7)
%       
%       fig3 = plotAmplitudeVectorFeild(-1, 4, -1, 1, 0.3, 1, 0.8, 3)
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Handle Arguments

if nargin < 4 || nargin > 8  
    disp('Error running plotAmplitudeVectorFeild')
    disp('Please provide the right number of arguments')
    disp('Function Call: ')
    disp('plotAmplitudeVectorFeild(alpha, beta1, beta2, epsilon)')
    disp('or')
    disp('plotAmplitudeVectorFeild(alpha, beta1, beta2, epsilon, F, f_osc, f_input)')
    fprintf('\n')
    disp('type "help plotAmplitudeVectorFeild" for more details')
    fprintf('\n')
    return
end

if nargin == 4
    f_osc = 1;
    f_input = 0.0;
    F = 0.0;
    figure_number = 1;
end

if nargin == 7
    figure_number = 1;
end

% make sure the input is properly zeroed out for autonomous osc
if F && ~f_input
    display('Both F and f_input must be zero for autonomous oscillators')
    display('F and f_input have been set to zero')
    f_input = 0;
    F = 0;
elseif  f_input && ~F
    display('Both F and f_input must be zero for autonomous oscillators')
    display('F and f_input have been set to zero')
    f_input = 0;
    F = 0;
end 

% Check input arguments
if F < 0
  error('F > 0')
end

if beta2 == 0 && epsilon ~= 0
    disp('set epsilon = 0 because beta2 = 0')
    epsilon = 0; 
end

if epsilon == 0 && beta2 ~= 0
    disp('set beta2 = 0 because epsilon = 0')
    beta2 = 0;
end

%% Get Param Regime and FPs

% Regime:
%   critical Hopf = 1
%   supercritical Hopf = 2
%   supercritical DLC = 3
%   subcritical DLC = 4
%   unknown = 0

% Stability Type:
%   stable node = 1 
%   stable spiral = 2
%   unstable node = 3
%   unstable spiral = 4
%   a saddle point = 5 
%   unknown = 0

[r_star, ~, stability_type, regime] = getFP(f_osc, f_input, alpha, beta1, beta2, epsilon, F, 0);

numFPs = length(r_star);
FPs = r_star;

%% sanity check

% Critical Hopf
if regime == 1 
    assert(numFPs == 1,'Critical Hopf has more than one FP')
% Supercritical Hopf
elseif regime == 2 % Supercritical Hopf
    assert(numFPs == 2,'Critical Hopf has more than one FP')
% Supercritical DLC
elseif regime == 3 
    assert(numFPs == 3,'Supercritical DLC does not have 3 FPs')
% Subcritical DLC
elseif regime == 4 
    assert(numFPs == 1,'Subcritical DLC has more than one FP')
% Unknown
elseif regime == 0 
    display('The oscillators parameter regime is unknown plot was not generated')
    return
end

%% Plot 

for i = 1:numFPs % loop through fixed points
    % plot
    plot_obj = figure(figure_number);
    if stability_type(i) == 1
        plot(FPs(i),0,'o','MarkerFaceColor', orange,'MarkerSize',10);
        hold on
    elseif stability_type(i) == 2
        plot(FPs(i),0,'o','MarkerFaceColor', 'y','MarkerSize',10);
        hold on
    elseif stability_type(i) == 3
        plot(FPs(i),0,'o','MarkerFaceColor', 'g','MarkerSize',10);
        hold on
    elseif stability_type(i) == 4
        plot(FPs(i),0,'o','MarkerFaceColor', 'c','MarkerSize',10);
        hold on
    elseif stability_type(i) == 5
        plot(FPs(i),0,'o','MarkerFaceColor', 'm','MarkerSize',10);
        hold on
    else
        display('FP could not be identified')
    end
end

function [C] = orange()
    C = [1 .5 0];
end

end