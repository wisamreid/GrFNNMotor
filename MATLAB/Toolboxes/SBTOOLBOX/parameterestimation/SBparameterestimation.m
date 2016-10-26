function [output] = SBparameterestimation(varargin)
% SBparameterestimation: performs estimation of selected parameters of the
% model
%
% USAGE:
% ======
% [output] = SBparameterestimation(model, parameters, data, options)
% [output] = SBparameterestimation(model, parameters, measurements, stimuli, options)
%
% model: SBmodel for which parameters are to be estimated
% parameters: structure with information about the parameters to estimate
%               parameters.names: cell array with names of parameters to estimate
%               parameters.initialValues: vector with initial guesses for parameters
%               parameters.signs: vector indicating limitations for the sign of given
%                   parameters (1 = pos, -1 = neg, 0 = can be both)
% data: SBdata object containing information about measurements and eventually also stimuli
% measurements: structure with measurement information
%               measurements.names: cell array with names of measured
%                   elements (states, variable, and/or reaction rates) to 
%                   use for estimation
%               measurements.values: matrix containing the measurement
%                   information. One row per time instant, one column per
%                   measured component. Values need to be defined for each
%                   time step in measurements.time (so far).
%               measurements.time: vector with time instants of
%                   measurements
% stimuli: structure with information about stimuli during the experiment
%               stimuli.names: cell array with names of stimulated
%                   elements (parameters) during experiment
%               stimuli.values: matrix containing the stimuli
%                   information. One row per time instant, one column per
%                   stimulated element. Values need to be defined for each
%                   time step in stimuli.time (so far).
%               stimuli.time: vector with time instants of stimuli
% options: structure with options for the parameter estimation
%               options.integrator: string with integrator to use for simulations
%               options.integratorOptions: options passed to integrator
%               options.initialConditions: initial conditions to start simulations from
%               options.optimizeInitialConditions: =0: initial conditions
%                      are not optimized, =1: initial conditions are optimized
%               options.plotResult: =0: do not plot final result, =1: plot
%                      final result
%               options.optimizer: string with optimization routine to use
%               options.optimizerOptions: options to pass to optimizer
%                   (please have a look at the functions "simplexSB" and
%                   "simannealingSB" to learn more about the possible
%                   options) You can also use other optimization routines,
%                   as long as the calling syntax is the same (e.g., fminsearch 
%                   from the MATLAB optimization toolbox).
%               options.optimizerCostFunction: string with name of cost function to
%                   use for estimation
%               options.optimizerOutputFunction: string with name of the output function that 
%                   is used to visualize the progress of the optimization
%               options.stimuliInterpolation: interpolation of stimuli values. 
%                   'linear' or 'nearest' are the options. Please have a
%                   look at the MATLAB function "interp1" to learn more
%                   about it.
%
% DEFAULT VALUES:
% ===============
% parameters.signs: assuming the signs of the parameters in the initial guess
% options.integrator: 'ode15s'
% options.integratorOptions: [] (default options)
% options.initialConditions: initial conditions stored in the model
% options.optimizeInitialConditions: 0 (no optimization of initial conditions)
% options.plotResult: 1 (plot result)
% options.optimizer: 'simplexSB'
% options.optimizerOptions: [] (default options)
% options.optimizerCostFunction: 'costFunctionParamEstSOS'
% options.optimizerOutputFunction: 'outputFunctionParamEst'
% options.stimuliInterpolation: 'linear'
%
% Output Arguments:
% =================
% output: structure with the results of the estimation
% output.success: 0=no success, 1=success
% output.parameterNames: names of the identified parameters (same as input
%       variable "parameters.names")
% output.parameterValues: vector with estimated parameter values
% output.initialconditionsOptimized: 'No' = initial conditions not
%       optimized, "Yes" = they have been optimized
% output.initialconditionsOptimal: optimized initial conditions if
%       optimized, otherwise the initial conditions are returned unchanged
% output.model: model in which the estimated parameters and optimized 
%       initial conditions have been set

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLOBAL VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global dataStructureCostFunction
global stimuliValues stimuliTime

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK IF SBMODEL OR FILENAME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp('SBmodel',class(varargin{1})),
    error('This function can not be applied to ODE files.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLE VARIABLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 4,
    model = varargin{1};
    parameters = varargin{2};
    data = varargin{3};
    options = varargin{4};
    [time, measurementNames, measurementValues, ...
        stimuliNames, stimuliValues] = SBgetdata(data);
    measurements = [];
    measurements.names = measurementNames;
    measurements.values = measurementValues;
    measurements.time = time;
    stimuli = [];
    stimuli.names = stimuliNames;
    stimuli.values = stimuliValues;
    stimuli.time = time;
    % check that all time points defined
    if ~isempty(find(isnan(measurementValues)==1)),
        error('Not all data points defined for all time points. Handling of that is not yet implemented.');
    end
    if ~isempty(find(isnan(stimuliValues)==1)),
        error('Not all data points defined for all time points. Handling of that is not yet implemented.');
    end
elseif nargin == 5,
    model = varargin{1};
    parameters = varargin{2};
    measurements = varargin{3};
    stimuli = varargin{4};
    options = varargin{5};
else
    error('Incorrect number of input arguments');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isfield(parameters,'names'), parameterNamesIdent = parameters.names; else error('No parameter names defined.'); end
if isfield(parameters,'initialValues'), parameterValuesStart = parameters.initialValues; else error('No parameter initial values defined.'); end
if isfield(parameters,'signs'), parameterSigns = parameters.signs; else parameterSigns = sign(parameterValuesStart); end
if isfield(measurements,'names'), measurementNames = measurements.names; else error('No measurement names defined.'); end
if isfield(measurements,'values'), measurementValues = measurements.values; else error('No measurement values defined.'); end
if isfield(measurements,'time'), measurementTime = measurements.time; else error('No measurement time vector defined.'); end
if isfield(stimuli,'names'), stimuliNames = stimuli.names; else stimuliNames = {}; end
if isfield(stimuli,'values'), stimuliValues = stimuli.values; else stimuliValues = []; end
if isfield(stimuli,'time'), stimuliTime = stimuli.time; else stimuliTime = []; end
if isfield(options,'integrator'), optionsIntegrator = options.integrator; else optionsIntegrator = 'ode15s'; end
if isfield(options,'integratorOptions'), optionsIntegratorOptions = options.integratorOptions; else optionsIntegratorOptions = []; end
if isfield(options,'initialConditions'), optionsInitialConditions = options.initialConditions; else optionsInitialConditions = []; end
if isfield(options,'optimizeInitialConditions'), optionsOptimizeInitialConditions = options.optimizeInitialConditions; else optionsOptimizeInitialConditions = 0; end
if isfield(options,'plotResult'), optionsPlotResult = options.plotResult; else optionsPlotResult = 1; end
if isfield(options,'optimizer'), optionsOptimizer = options.optimizer; else optionsOptimizer = 'simplexSB'; end
if isfield(options,'optimizerOptions'), optionsOptimizerOptions = options.optimizerOptions; else optionsOptimizerOptions = []; end
if isfield(options,'optimizerCostFunction'), optionsOptimizerCostFunction = options.optimizerCostFunction; else optionsOptimizerCostFunction = 'costFunctionParamEstSOS'; end
if isfield(options,'optimizerOutputFunction'), optionsOptimizerOutputFunction = options.optimizerOutputFunction; else optionsOptimizerOutputFunction = 'outputFunctionParamEst'; end
if isfield(options,'stimuliInterpolation'), stimuliInterpolation = options.stimuliInterpolation; else stimuliInterpolation = ''; end

% if single identifiers given as strings then convert them to cell arrays
% (parameter names, measurement names, and stimuli names)
if strcmp(class(parameterNamesIdent),'char'),
    parameterNamesIdent = {parameterNamesIdent};
end
if strcmp(class(measurementNames),'char'),
    measurementNames = {measurementNames};
end
if strcmp(class(stimuliNames),'char'),
    stimuliNames = {stimuliNames};
end
% if no initial parameter values are given then use the 
% ones stored in the model
if isempty(parameterValuesStart),
    [names, values] = SBparameters(model);
    for k1 = 1:length(parameterNamesIdent),
        for k2 = 1:length(names),
            if strcmp(names{k2},parameterNamesIdent{k1}),
                parameterValuesStart(k1) = values(k2);
                break;
            end
        end           
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK INPUT ERRORS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% parameters names for parameters to identify need to be given
if length(parameterNamesIdent) == 0,
    error('No parameters names given.');
end
% check if same number of parameter names and initial 
% parameter values given
if length(parameterValuesStart) ~= length(parameterNamesIdent),
    error('Inconsistent number of parameter names and initial parameter values.');
end
% if parameterSigns empty then assume parameters should keep their signs
if isempty(parameterSigns),
    parameterSigns = sign(parameterValuesStart);
end
% check parameter signs setting
if length(parameterSigns) ~= length(parameterNamesIdent),
    error('Inconsistent number of parameter names and parameter signs setting.');
end
% check that sign settings correspond to initial parameter signs
test = parameterSigns(:).*parameterValuesStart(:);
if ~isempty(find(test < 0)), 
    error('Conflict between initial parameter values and parameter sign settings');
end
% if initial parameter is zero but sign setting -1 or +1 then just warn
for k = 1:length(parameterValuesStart),
    if parameterValuesStart(k) == 0,
        disp(sprintf('In case the sign of a parameter is of importance you should avoid\nto set initial values of parameters to 0.\nYou might want to consider a different starting guess.'));
    end
end
% measurement names need to be given and the number
% should be the same as measurement data available
% also the measurements should have the same number of time
% steps as the measurement time vector
if length(measurementNames) == 0,
    error('No measurement names given.');
end
if length(measurementNames) > 1,
    if length(measurementNames) ~= size(measurementValues,2),
        error('Inconsistent number of measurement names and measurement data.');
    end
    if length(measurementTime) ~= size(measurementValues,1),
        error('Inconsistent number of measurement time steps and measurement data.');
    end
elseif length(measurementNames) == 1,
    if min(size(measurementValues)) ~= 1,
        error('Inconsistent number of measurement names and measurement data.');
    end
    if max(size(measurementTime)) ~= max(size(measurementValues)),
        error('Inconsistent number of measurement time steps and measurement data.');
    end    
end
% stimuli parameter names need NOT to be given. the number of given
% stimuli should be the same as stimuli data available
% also the stimuli should have the same number of time
% steps as the stimuli time vector
if length(stimuliNames) > 1,
    if length(stimuliNames) ~= size(stimuliValues,2),
        error('Inconsistent number of stimuli names and stimuli data.');
    end
    if length(stimuliTime) ~= size(stimuliValues,1),
        error('Inconsistent number of stimuli time steps and stimuli data.');
    end
elseif length(stimuliNames) == 1,
    if min(size(stimuliValues)) ~= 1,
        error('Inconsistent number of stimuli names and stimuli data.');
    end
    if max(size(stimuliTime)) ~= max(size(stimuliValues)),
        error('Inconsistent number of stimuli time steps and stimuli data.');
    end    
    % if just one stimuli then make sure the values and time are given as
    % column vectors
    stimuliTime = stimuliTime(:);
    stimuliValues = stimuliValues(:);
end
% check if initial conditions given and if the length of the vector matches
% the number of states in the model
if length(optionsInitialConditions) > 0 && length(optionsInitialConditions) ~= length(SBinitialconditions(model)),
    error('Number of initial conditions does not match number of states in the model.');
end
% check stimuli interpolation setting
if length(stimuliNames) >= 1,
    if ~strcmp(stimuliInterpolation, 'linear') && ~strcmp(stimuliInterpolation, 'nearest'),
        stimuliInterpolation = 'linear';
        disp('No interpolation method for stimuli specified. Assuming linear interpolation.');
        disp(' ');
    end
end
% check that the time vector for the stimuli starts not later and ends not 
% earlier than the time vector of the measurements!
if ~isempty(stimuliTime),
    warning = 0;
    if stimuliTime(1) > measurementTime(1),
        stimuliTime = [measurementTime(1); stimuliTime];
        stimuliValues = [stimuliValues(1,:); stimuliValues];
        warning = 1;
    end
    if stimuliTime(end) < measurementTime(end),
        stimuliTime = [stimuliTime; measurementTime(end)];
        stimuliValues = [stimuliValues; stimuliValues(end,:)];
        warning = 1;
    end
    if warning == 1,
        disp(sprintf('Warning: stimuli data not defined over whole time range of\nmeasurement data. The function will extrapolate as if the stimuli\nare constant until the first data point and from the last data point\nto the end time of the measurements.\n'));
        disp('Press a key to continue');
        pause;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK IF GIVEN PARAMETER NAMES EXIST IN MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (parameters)
parametersModel = SBparameters(model);
errorMsg = '';
for k1 = 1:length(parameterNamesIdent),
    found = 0;
    for k2 = 1:length(parametersModel),
        if strcmp(parametersModel{k2},parameterNamesIdent{k1}),
            found = 1;
            break;
        end
    end
    if found == 0,
        errorMsg = sprintf('%sParameter %s is not present in the model\n',errorMsg,parameterNamesIdent{k1});
    end
end
% error message in case that at least one parameter not present in the model
if length(errorMsg) > 0,
    error(errorMsg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK IF MEASUREMENT ELEMENTS PRESENT IN MODEL (ALSO GET THEIR INDICES)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (states, variables, reaction rates)
errorMsg = '';
statesModel = SBstates(model);
variablesModel = SBvariables(model);
reactionsModel = SBreactions(model);
% indices are stored as: first row: index of element in measurement
%                        second row: index of element in element vector
stateIndices = [];
variableIndices = [];
reactionIndices = [];
for k1 = 1:length(measurementNames),
    found = 0;
    for k2 = 1:length(statesModel),
        if strcmp(statesModel{k2},measurementNames{k1}),
            found = 1;
            stateIndices = [stateIndices [k1; k2]];
            break;
        end
    end
    if found == 0,
        for k2 = 1:length(variablesModel),
            if strcmp(variablesModel{k2},measurementNames{k1}),
                found = 1;
                variableIndices = [variableIndices [k1; k2]];
                break;
            end
        end
    end
    if found == 0,
        for k2 = 1:length(reactionsModel),
            if strcmp(reactionsModel{k2},measurementNames{k1}),
                found = 1;
                reactionIndices = [reactionIndices [k1; k2]];
                break;
            end
        end
    end
    if found == 0,
        errorMsg = sprintf('%sMeasurement %s is not present in the model\n',errorMsg,measurementNames{k1});
    end
end
% error message in case that at least one measurement not present in the model
if length(errorMsg) > 0,
    error(errorMsg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK IF STIMULI PARAMETERS PRESENT IN MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (parameters)
parametersModel = SBparameters(model);
errorMsg = '';
for k1 = 1:length(stimuliNames),
    found = 0;
    for k2 = 1:length(parametersModel),
        if strcmp(parametersModel{k2},stimuliNames{k1}),
            found = 1;
            break;
        end
    end
    if found == 0,
        errorMsg = sprintf('%sStimulus parameter %s is not present in the model\n',errorMsg,stimuliNames{k1});
    end
end
% error message in case that at least one stimulus parameter not present in the model
if length(errorMsg) > 0,
    error(errorMsg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MODIFY THE MODEL FOR STIMULI USE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the parameters used for stimuli are deleted from the model and replaced
% by functions that do interpolation between the stimuli values
% first delete the parameters
modelIdent = deleteparametersSB(model,stimuliNames);
% then replace the deleted parameters in the equations by functions
for k = 1:length(stimuliNames),
    modelIdent = replaceelementSB(modelIdent,stimuliNames{k},strcat(stimuliNames{k},'(time)'));
end
% construct the matlab functions for the stimuli parameters
matlabFunctions = '';
for k = 1:length(stimuliNames),
    matlabFunctions = sprintf('%sfunction [%s] = %s(time)\n',matlabFunctions,stimuliNames{k},stimuliNames{k});
    matlabFunctions = sprintf('%sglobal stimuliValues stimuliTime\n',matlabFunctions);
    if strcmp(stimuliInterpolation, 'linear'),
        matlabFunctions = sprintf('%s%s = interp1(stimuliTime, stimuliValues(:,%d)'', time,''linear'');\n',matlabFunctions,stimuliNames{k},k);
    else
        matlabFunctions = sprintf('%s%s = interp1(stimuliTime, stimuliValues(:,%d)'', time,''nearest'');\n',matlabFunctions,stimuliNames{k},k);
    end
    matlabFunctions = sprintf('%sreturn\n\n',matlabFunctions);
end
% add matlabFunctions to model
modelStructure = SBstruct(modelIdent);
oldfunctions = modelStructure.functionsMATLAB;
newfunctions = sprintf('%s\n%s',oldfunctions, matlabFunctions);
modelStructure.functionsMATLAB = newfunctions;
modelIdent = SBmodel(modelStructure);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REORDER THE MEASUREMENTS AND CONSTRUCT INDEX VECTORS FOR SIMULATION DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(stateIndices)
    stateIndices = sortrows(stateIndices',2)';
    stateMeasurements = measurementValues(:,stateIndices(1,:));
    stateIndices = stateIndices(2,:);
else
    stateMeasurements = [];
end
if ~isempty(variableIndices)
    variableIndices = sortrows(variableIndices',2)';
    variableMeasurements = measurementValues(:,variableIndices(1,:));
    variableIndices = variableIndices(2,:);
else
    variableMeasurements = [];
end
if ~isempty(reactionIndices)
    reactionIndices = sortrows(reactionIndices',2)';
    reactionMeasurements = measurementValues(:,reactionIndices(1,:));
    reactionIndices = reactionIndices(2,:);
else
    reactionMeasurements = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCT DATA STRUCTURE FOR COST FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataStructureCostFunction = [];
dataStructureCostFunction.model = modelIdent;
dataStructureCostFunction.time = measurementTime;
if isempty(optionsInitialConditions),
    dataStructureCostFunction.initialConditions = SBinitialconditions(modelIdent);
else
    dataStructureCostFunction.initialConditions = optionsInitialConditions;
end
dataStructureCostFunction.optimizeInitialConditions = optionsOptimizeInitialConditions;
dataStructureCostFunction.parameters = parameterNamesIdent;
dataStructureCostFunction.integrator = optionsIntegrator;
dataStructureCostFunction.integratorOptions = optionsIntegratorOptions;
dataStructureCostFunction.stateMeasurements = stateMeasurements;
dataStructureCostFunction.stateIndices = stateIndices;
dataStructureCostFunction.variableMeasurements = variableMeasurements;
dataStructureCostFunction.variableIndices = variableIndices;
dataStructureCostFunction.reactionMeasurements = reactionMeasurements;
dataStructureCostFunction.reactionIndices = reactionIndices;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLE THE CASE OF OPTIMIZATION OF INITIAL CONDITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if optionsOptimizeInitialConditions == 1,
    startValues = [dataStructureCostFunction.initialConditions(:)' parameterValuesStart(:)'];
else
    startValues = parameterValuesStart(:)';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK OPTIMIZER OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if no parameter signs given then use the signs of the initial guesses for
% the parameters. 
if strcmp(optionsOptimizer,'simplexSB'),
    if length(optionsOptimizerOptions) < 6,
        optionsOptimizerOptions{6} = sign(parameterValuesStart(:)');
    elseif isempty(optionsOptimizerOptions{6}),
        optionsOptimizerOptions{6} = sign(parameterValuesStart(:)');
    end
    optionsOptimizerOptions{7} = optionsOptimizerOutputFunction;
end
if strcmp(optionsOptimizer,'simannealingSB'),
    if length(optionsOptimizerOptions) < 10,
        optionsOptimizerOptions{10} = sign(parameterValuesStart(:)');
    elseif isempty(optionsOptimizerOptions{10}),
        optionsOptimizerOptions{10} = sign(parameterValuesStart(:)');
    end
    optionsOptimizerOptions{11} = optionsOptimizerOutputFunction;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO THE OPTIMIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Xoptimal,FVAL,EXITFLAG]= feval(optionsOptimizer,optionsOptimizerCostFunction,startValues,optionsOptimizerOptions);

if dataStructureCostFunction.optimizeInitialConditions == 0,
    initialConditionsOptimal = dataStructureCostFunction.initialConditions;
    parametersOptimal = Xoptimal;
else
    initialConditionsOptimal = Xoptimal(1:length(dataStructureCostFunction.initialConditions));
    parametersOptimal = Xoptimal(length(dataStructureCostFunction.initialConditions)+1:end);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS RETURN PARAMETERS AND DISPLAY ...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelOptimal = SBparameters(modelIdent,parameterNamesIdent,parametersOptimal);
modelOptimal = SBinitialconditions(modelOptimal, initialConditionsOptimal);
if EXITFLAG == 1,
    success = 'Yes';
else
    success = 'No';
end
output.success = success;
output.parameterNames = parameterNamesIdent;
output.parameterValues = parametersOptimal;
if dataStructureCostFunction.optimizeInitialConditions == 0,
    output.initialconditionsOptimized = 'No';
else
    output.initialconditionsOptimized = 'Yes';
end
output.initialconditionsOptimal = initialConditionsOptimal;
output.model = modelOptimal;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY PLOT ON REQUEST
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if optionsPlotResult == 1,
    simulationData = SBsimulate(modelOptimal, ...
        dataStructureCostFunction.integrator, dataStructureCostFunction.time, ...
        initialConditionsOptimal, dataStructureCostFunction.integratorOptions);
    % get the simulated values
    stateSimulations = simulationData.statevalues(:,dataStructureCostFunction.stateIndices);
    variableSimulations = simulationData.variablevalues(:,dataStructureCostFunction.variableIndices);
    reactionSimulations = simulationData.reactionvalues(:,dataStructureCostFunction.reactionIndices);
    % plot results
    figure(1); clf;
    simulationResults = real([stateSimulations variableSimulations reactionSimulations]);
    measurementCompare = [stateMeasurements variableMeasurements reactionMeasurements];
    plot(measurementTime, simulationResults,'-'); hold on;
    plot(measurementTime, measurementCompare,'--o'); hold on;
end
return

		


