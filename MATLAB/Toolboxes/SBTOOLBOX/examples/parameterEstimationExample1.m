clc; clear all; close all;
echo on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER ESTIMATION EXAMPLE 1
% Note: Please take this file as an example on how to use the different
% functions, not as a real life example.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define the model to use (Novak-Tyson Cell Cycle model)
model = SBmodel('parameterEstimationExampleModel.txt');
% Press a key
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PERFORMING IN SILICO EXPERIMENTS TO COLLECT MEASUREMENT DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the measurements (here all states are measured)
% The sampling time is chosen to be 5min, measurements start at time 0 and
% end at time 100.
measurements = {'states', [5 100]};
% Press a key
pause;

% Define the noise added to the measurements {offset,variance}
% values are in percent of the mean values of the measurements 
% during the simulation time
noise = {0,0.1}; 
% Press a key
pause;

% No perturbations are applied to the system during the experiment
perturbation = {};
% Press a key
pause;

% Running the insilico experiment
[datanoise, data, noiseoffset, noisevariance, experimentmodel] = SBexperiment(model, measurements, perturbation, noise);
% Press a key
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VISUALUIZE MEASUREMENT DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SBvisualizedata(datanoise);
% Press a key
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETER ESTIMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the names of the parameters to estimate
parameters.names = {'Ka','Kb','Kc'};
% Press a key
pause;

% Define the initial values (the true values are: [0.1 1 0.01]
parameters.initialValues = [10  100 20];
% Press a key
pause;

% The parameters need to have positive values
parameters.signs = [1 1 1];
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
result = SBparameterestimation(model, parameters, datanoise, options)
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
