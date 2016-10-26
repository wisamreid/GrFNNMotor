clc; clear all;
echo on
% Example script illustrating the use of the "SBnetworkident" function for
% for the identification of a local connectivity matrix for measured
% components.
%
% Due to lack of real experimental data we will use here in silico
% generated measurement data. This data is generated using the SBexperiment
% function. The model used is a 4 gene network.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEFINING THE MODEL TO USE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
model = SBmodel('kholodenko4gene.txt');
% Determine the steady-state of the system
[ss,res] = SBsteadystate(model);
% Update the model with steady-state information (in order to start the
% experiment while system is in its steady-state)
model = SBinitialconditions(model,ss);
% Press a key
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PERFORMING IN SILICO EXPERIMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define the measurements (here all states are measured)
% The sampling time is chosen to be 0.01, measurements start at time 0 and
% end at time 0.2.
measurements = {'states', [0.01 0.2]};
% Press a key
pause;

% Define the noise added to the measurements {offset,variance}
% values are in percent of the mean values of the measurements 
% during the simulation time
noise = {0,0}; 
% Press a key
pause;

% Define the perturbation to be applied to the system
perturbationSize = 50; % (in percent)
% Experiment 1) Perturbation for first experiment [time instant, perturbation, absRel]
perturbation1 = {'Vs1','step',[0.1 perturbationSize 1]};
% Experiment 2) Perturbation for second experiment [time instant, perturbation, absRel]
perturbation2 = {'Vs2','step',[0.1 perturbationSize 1]};
% Experiment 3) Perturbation for third experiment [time instant, perturbation, absRel]
perturbation3 = {'Vs3','step',[0.1 perturbationSize 1]};
% Experiment 4) Perturbation for fourth experiment [time instant, perturbation, absRel]
perturbation4 = {'Vs4','step',[0.1 perturbationSize 1]};
% Press a key
pause;

% running the 4 experiments
[datanoise1, data1, noiseoffset, noisevariance, experimentmodel] = SBexperiment(model, measurements, perturbation1, noise);
[datanoise2, data2, noiseoffset, noisevariance, experimentmodel] = SBexperiment(model, measurements, perturbation2, noise);
[datanoise3, data3, noiseoffset, noisevariance, experimentmodel] = SBexperiment(model, measurements, perturbation3, noise);
[datanoise4, data4, noiseoffset, noisevariance, experimentmodel] = SBexperiment(model, measurements, perturbation4, noise);
% Press a key
pause;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PERFORMING THE NETWORK IDENTIFICATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EstimatedConnectivityMatrix = SBnetworkident(datanoise1, datanoise2, datanoise3, datanoise4, [12 12 12 12],[5 5 5 5])
% Press a key
pause;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARE THE RESULT TO THE JACOBIAN OF THE SYSTEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Under the following assumptions the network connectivity matrix N
% corresponds to an estimation of the Jacobian of the system at the
% considered steady-state: 
%   - sampling time chosen small enough
%   - all components in the network are measured
%   - perturbations chosen small enough such that the response of the
%     system can be seen as linear
Jacobian = SBjacobian(model,SBsteadystate(model))
% Press a key
pause;

% We can see that even if the perturbation size has been chosen relatively
% large (50% activation), the estimated network connectivity matrix is
% reasonably close to the Jacobian. Choosing smaller perturbations this
% result gets of course better.
echo off