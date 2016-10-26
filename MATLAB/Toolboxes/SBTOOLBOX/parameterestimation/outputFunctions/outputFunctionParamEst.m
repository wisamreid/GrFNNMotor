function [ output_args ] = outputFunctionParamEst(X,FVAL)
% outputFunctionParamEst: Default output function for displaying
% intermediate results during parameter estimation. The global variable
% "dataStructureCostFunction" is defined in the SBparameterestimation
% function and can be read out and used here.

% Information:
% ============
% Systems Biology Toolbox for MATLAB
% Copyright (C) 2005-2007 Fraunhofer-Chalmers Research Centre
% Author of the toolbox: Henning Schmidt, henning@hschmidt.de
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details. 
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
% USA.

global dataStructureCostFunction
% handle the case when initial conditions are to be optimized additionally
if dataStructureCostFunction.optimizeInitialConditions == 0,
    initialConditions = dataStructureCostFunction.initialConditions;
    parameters = X;
else
    initialConditions = X(1:length(dataStructureCostFunction.initialConditions));
    parameters = X(length(dataStructureCostFunction.initialConditions)+1:end);
end
% set the given parameters in the model
model = SBparameters(dataStructureCostFunction.model, dataStructureCostFunction.parameters, parameters);
% simulate the model
simulationData = SBsimulate(model, ...
    dataStructureCostFunction.integrator, dataStructureCostFunction.time, ...
    initialConditions, dataStructureCostFunction.integratorOptions);
% get the simulated values for display
stateSimulations = simulationData.statevalues(:,dataStructureCostFunction.stateIndices);
variableSimulations = simulationData.variablevalues(:,dataStructureCostFunction.variableIndices);
reactionSimulations = simulationData.reactionvalues(:,dataStructureCostFunction.reactionIndices);
elementSimulations = real([stateSimulations, variableSimulations, reactionSimulations]);
% get the names of measured and simulated elements
stateNames = SBstates(dataStructureCostFunction.model);
variableNames = SBvariables(dataStructureCostFunction.model);
reactionNames = SBreactions(dataStructureCostFunction.model);
stateNames = stateNames(dataStructureCostFunction.stateIndices);
variableNames = variableNames(dataStructureCostFunction.variableIndices);
reactionNames = reactionNames(dataStructureCostFunction.reactionIndices);
elementNames = {stateNames{:} variableNames{:} reactionNames{:}};
% construct measurement element
elementMeasurements = [dataStructureCostFunction.stateMeasurements, dataStructureCostFunction.variableMeasurements, dataStructureCostFunction.reactionMeasurements];
% plot results
figure(1); clf;
plot(dataStructureCostFunction.time, elementSimulations); hold on;
plot(dataStructureCostFunction.time, elementMeasurements,'--o');
legend(elementNames);
drawnow;
return