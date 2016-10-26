function [varargout] = costFunctionParamEstSOS(X)
% costFunctionParamEstSOS: cost function for use for parameter estimation.
% The cost is determined by the sum of squares of the differences of
% simulated and measurement data. The function allows to determine the cost
% both for changed parameters and initial conditions.
% The global variable "dataStructureCostFunction" is initialized by the 
% SBparameterestimation function.
%
% The function has a variable number of output arguments:
%   output = [cost]
%   output = [cost, constraints]
% The variable "cost" should be a scalar value defining the cost of a
% certain set of parameters.
% The variable "constraints" can be used to pass information about
% constraints to the optimization algorithm. Not all optimizers are able to
% handle constraints / expect the corresponding output argument -  and thus
% this variable is optional. 
%
% The user has to add the equations that are needed for "constraints". In
% this default costfunction, if called with two output arguments,
% constraints = [].

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cost
%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
% get the simulated values to compare to the measurements
stateSimulations = real(simulationData.statevalues(:,dataStructureCostFunction.stateIndices));
variableSimulations = real(simulationData.variablevalues(:,dataStructureCostFunction.variableIndices));
reactionSimulations = real(simulationData.reactionvalues(:,dataStructureCostFunction.reactionIndices));
if isempty(stateSimulations), stateSimulations = []; end
if isempty(variableSimulations), variableSimulations = []; end
if isempty(reactionSimulations), reactionSimulations = []; end
% determine differences between measurement and simulation data
stateDifference = stateSimulations - dataStructureCostFunction.stateMeasurements;
variableDifference = variableSimulations - dataStructureCostFunction.variableMeasurements;
reactionDifference = reactionSimulations - dataStructureCostFunction.reactionMeasurements;
% determine cost by summing the squared differences
stateCost = sum(sum(stateDifference.^2));
variableCost = sum(sum(variableDifference.^2));
reactionCost = sum(sum(reactionDifference.^2));
% determine total cost by adding up the three
cost = stateCost + variableCost + reactionCost;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constraints
%%%%%%%%%%%%%%%%%%%%%%%%%%%
constraints = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout == 1,
    varargout{1} = cost;
elseif nargout == 2,
    varargout{1} = cost;
    varargout{2} = constraints;
end
return