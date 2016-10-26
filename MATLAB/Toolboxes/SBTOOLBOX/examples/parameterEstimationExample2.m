clc; clear all; close all;
echo on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER ESTIMATION EXAMPLE 2
% Note: Please take this file as an example on how to use the different
% functions, not as a real life example.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define the model to use (Novak-Tyson Cell Cycle model)
model = SBmodel('parameterEstimationExampleModel.txt');
% Press a key
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATION OF MEASUREMENT DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the measurements (here all states are measured)
% The sampling time is chosen to be 5min, measurements start at time 0 and
% end at time 100.
time = [0:5:100];
% Press a key
pause;

% Simulation of the model to obtain measurement data
measurementData = SBsimulate(model,'ode23s',time,SBinitialconditions(model));
% Press a key
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER ESTIMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the names of the parameters to estimate
parameters.names = {'Ka','Kb','Kc'};
% Define the initial values (the true values are: [0.1 1 0.01]
parameters.initialValues = [10  100 20];
% The parameters need to have positive values
parameters.signs = [1 1 1];
% Press a key
pause;

% Use all state measurements for the parameter estimation
measurements.names = SBstates(model);  % measure all states
% Get the measurement values and add some random noise
noise = 0.1;
measurements.values = measurementData.statevalues.*(1+noise*randn(size(measurementData.statevalues)));
% Time vector for measurements is the same as the time vector for the
% simulation (insilico experiment) above
measurements.time = time;
% Press a key
pause;

% No stimuli used for perturbation during the experiment
stimuli = [];
% Press a key
pause;

% Use simulated annealing as optimization method to find
% global optimum
options.optimizer = 'simannealingSB'; 
% Press a key
pause;

% Use default optimizer options
options.optimizerOptions = {};
% Press a key
pause;

% Start the parameter estimation
result = SBparameterestimation(model, parameters, measurements, stimuli, options)
% Press a key
pause;

% The algorithm should in many cases be able to converge to parameter 
% values very close to the real values [0.1 1 0.01]. Due to the added
% noise, the estimate can not be perfect. Setting the added noise to 0 will
% lead to a perfect estimation of the parameters. In some cases the
% algorithm might get stuck in a local minimum.
result.parameterValues
% Press a key
pause;

% You can try a local optimization method by exchanging 'simannealingSB'
% against 'simplexSB' in the above code. The result in this example should
% always be a non optimal estimate for the parameters.




