function [datanoise,data,noiseoffset,noisevariance,experimentmodel] = SBexperiment(varargin)
% SBexperiment: This function allows to do in silico experiments. You can
% specify the things you want to measure (states, variables, and/or
% reaction rates), perturbation you want to apply to the system (in terms
% of parameter perturbations), and specify noise that is added to the
% measurements.
%
% The experiment starts always at the initial state that is specified in
% the model. If you want to start at a specific state you will have to set
% this state correclty in the model by using SBedit or the
% SBinitialconditions function.
%
% USAGE:
% ======
% [datanoise,data,noiseoffset,noisevariance,experimentmodel] = SBexperiment(model,measurements,perturbations)
% [datanoise,data,noiseoffset,noisevariance,experimentmodel] = SBexperiment(model,measurements,perturbations,noise)
% [datanoise,data,noiseoffset,noisevariance,experimentmodel] = SBexperiment(model,measurements,perturbations,noise,integratorMethod)
% [datanoise,data,noiseoffset,noisevariance,experimentmodel] = SBexperiment(model,measurements,perturbations,noise,integratorMethod,integratorOptions)
%
% model: SBmodel to perform the in silico experiment on
% measurements: cell-array defining the in-silico measurements that are to 
%   be done. The format for this cell-array is
%   
%       measurements = {name1, time1, name2, time2, ..., namen, timen}
%
%       name1-n: string with the name of an element in the model to be
%                measured. An element can be a state, a variable, or a
%                reaction rate. Instead of specifying a string with one
%                name it is also possible to specify name1-n as a
%                cell-array containing strings with element names. 
%                
%                Special cases for names: 
%                   'all': measuring all elements
%                   'states': measuring all states
%                   'variables': measuring all variables
%                   'reactions': measuring all reaction rates
%
%                The time information specified in the case of the special
%                names applies to all elements that are selected by the
%                special identifier. The same holds if several names are
%                specified in one name argument by using a cell-array.
%       
%       time1-n: vector with time information at which instants to take the 
%                measurement samples for the corresponding element.
%                If this vector contains 2 elements, the format is the
%                following: 
%                           timex = [deltaT endTime]
%
%                where deltaT is the sampling time and endTime the time 
%                until which to collect the measurements. If the time
%                vector contains more than 2 elements then these are
%                interpreted as the time steps at which to measure the
%                corresponding element.
%                
%                The time vector can be chosen differently for each element
%                to be measured.
%      
% perturbations: cell-array defining the perturbations to be applied to the
%   model during the experiment. The format of the cell-array is the
%   following:
%       
%       perturbations = {parameterName1, type1, values1, ..., parameterNamen, typen, valuesn}
%
%       parameterName1: string defining the name of the parameter to be
%                       changed by the perturbation.
%
%       type1-n: string defining the type of the perturbation. Depending on
%                the type of the perturbation the values1-n vectors contain
%                different data. The different types and the correponding
%                syntax of the values are:
%                
%                'step': performing a step in given parameter at a given
%                        time, defined by 'timeInstant'. The size of this
%                        step is defined by 'sizePerturbation'. Per default 
%                        this perturbation is assumed to be relative and
%                        that 'sizePerturbation' determines the
%                        perturbation in percent. Specifying 'absRel' the
%                        user can select if the perturbation is absolute
%                        or relative. absRel = 0: absolute, = 1: relative
%                        (percent) (default).
%
%                        valuesx = [timeInstant, sizePerturbation] 
%                        valuesx = [timeInstant, sizePerturbation, absRel] 
%
%                'pulse': performing a pulse in given parameter, changing
%                         its nominal value to a given value at a given
%                         time 'timInstantChange' and changing back to its
%                         nominal value at a given time 'timeInstantChangeBack'.
%                         The size of this step is defined by
%                         'sizePerturbation'. Per default this perturbation
%                         is assumed to be relative and that
%                         'sizePerturbation' determines the perturbation in
%                         percent. Specifying 'absRel' the user can select
%                         if the perturbation is absolute or relative.
%                         absRel = 0: absolute, = 1: relative percent) (default).
%
%                         valuesx = [timeInstantChange, timeInstantChangeBack, sizePerturbation]
%                         valuesx = [timeInstantChange, timeInstantChangeBack, sizePerturbation, absRel]
%
% noise:      cell-array defining the noise (normally distributed) to be
%             added to the measurements. The format for this input argument is:    
%   
%             noise = {offset, variance}
%             noise = {offset, variance, absRel}
%
%             offset:   either a scalar value determining the noise offset for
%                       all measurements, or a vector of the same length as
%                       the number of measurements allowing to specify a
%                       different noise offset for each measurement (same
%                       order as measurements)
%             variance: either a scalar value determining the noise variance for
%                       all measurements, or a vector of the same length as
%                       the number of measurements allowing to specify a
%                       different noise variance for each measurement (same
%                       order as measurements)
%             absRel:   =0: noise data is specified in absolute values
%                       =1: noise data is specified in relative values,
%                       that is here, in percent of the mean of the
%                       absolute values of the measurement values over the 
%                       entire experiment. (default setting) 
%
% DEFAULT VALUES:
% ===============
% noise: if not specified, no noise will be added to the measurements
% integratorMethod: 'ode23s'
% integratorOptions: [] (MATLAB default)
%
% Output Arguments:
% =================
% If no output arguments are specified the experiment results are
% displayed using the SBvisualizedata function. (The datanoise SBdata
% object is displayed, not that data SBobject)
%
% datanoise: SBdata object containing the in silico measurements (with
%   noise - in the case that no noise information is given, it will contain
%   the same data as the 'data' output variable).
% data: SBdata object containing the in silico measurements (without noise)
% noiseoffset: the noise offset used for generating the noise in absolut
%   values (not in relative values this means)
% noisevariance: the noise variance used for generating the noise in absolut
%   values (not in relative values this means)
% experimentmodel: SBmodel that has been used to generate the measurements

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
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK IF SBmodel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp('SBmodel',class(varargin{1})),
    error('Function only defined for SBmodels.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY INFORMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Performing in silico experiment');
notesText = 'In silico experiment';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS VARIABLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 3,
    model = varargin{1};
    measurements = varargin{2};
    perturbations = varargin{3};
    noise = {0,0};
    integratorMethod = 'ode23s';
    integratorOptions = [];
elseif nargin == 4,
    model = varargin{1};
    measurements = varargin{2};
    perturbations = varargin{3};
    noise = varargin{4};
    integratorMethod = 'ode23s';
    integratorOptions = [];
elseif nargin == 5,
    model = varargin{1};
    measurements = varargin{2};
    perturbations = varargin{3};
    noise = varargin{4};
    integratorMethod = varargin{5};
    integratorOptions = [];
elseif nargin == 6,
    model = varargin{1};
    measurements = varargin{2};
    perturbations = varargin{3};
    noise = varargin{4};
    integratorMethod = varargin{5};
    integratorOptions = varargin{6};
else 
    error('Incorrect number of input arguments.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS THE MEASUREMENT INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% cycle through all elements in the measurement vector and construct a list
% of the elements to be measured and the corresponding time vectors.
% Perform checks on the time vectors that they are acceptably defined.
% Construct a master time vector containing all the time instants that are
% specified for any measurement (this is then used as time vector for the
% simulation).
elementsMeasured = {};  % containing information about elements to be measured and corresponding time instants
masterTimeVector = [];
% check that even number of elements
if mod(length(measurements),2) ~= 0,
    error('Incorrect number of elements in the measurements input argument.');
end
elementsMeasuredIndex = 1;
for k = 1:length(measurements)/2,
    nameMeasurement = strtrim(measurements{2*k-1});
    timeMeasurement = measurements{2*k};
    % check time vector and expand it if needed
    if length(timeMeasurement) < 2,
        error('Time vector for measurement ''%s'' not correctly defined.',nameMeasurement);
    end
    if length(timeMeasurement) == 2,
        % check that sampling tim larger than zero
        if timeMeasurement(1) == 0,
            error('Sampling time for measurement ''%s'' is set to zero.',nameMeasurement);
        end
        timeMeasurement = [0:timeMeasurement(1):timeMeasurement(2)];
    end
    % check that the timevector is monotonically increasing
    test = find(timeMeasurement(2:end)-timeMeasurement(1:end-1) < 0);
    if ~isempty(test),
        error('Time vector for measurement ''%s'' is not monotonically increasing.',nameMeasurement);
    end
    % check that only positive time instants defined
    test = find(timeMeasurement < 0);
    if ~isempty(test),
        error('Negative time instants present in time vector for measurement ''%s''.',nameMeasurement);
    end
    % get the element names in the special cases
    specialCase = 0;
    if strcmp(nameMeasurement,'all'),
        specialCase = 1;
        specialNames1 = SBstates(model);
        specialNames2 = SBvariables(model);
        specialNames3 = SBreactions(model);
        specialNames = {specialNames1{:},specialNames2{:},specialNames3{:}};
    elseif strcmp(lower(nameMeasurement),'states'),
        specialCase = 1;
        specialNames = SBstates(model);
    elseif strcmp(lower(nameMeasurement),'variables'),
        specialCase = 1;
        specialNames = SBvariables(model);
    elseif strcmp(lower(nameMeasurement),'reactions'),
        specialCase = 1;
        specialNames = SBreactions(model);
    elseif strcmp(class(nameMeasurement),'cell'),
        % elements to measure given as cell-array of strings
        specialCase = 1;
        specialNames = nameMeasurement;
    end
    % add the time and name data to structure
    if specialCase == 0,
        elementsMeasured(elementsMeasuredIndex).name = nameMeasurement;
        elementsMeasured(elementsMeasuredIndex).time = timeMeasurement;
        masterTimeVector = unique([masterTimeVector(:); timeMeasurement(:)]);
        elementsMeasuredIndex = elementsMeasuredIndex + 1;
    else
        for k2 = 1:length(specialNames),
            elementsMeasured(elementsMeasuredIndex).name = specialNames{k2};
            elementsMeasured(elementsMeasuredIndex).time = timeMeasurement;
            masterTimeVector = unique([masterTimeVector(:); timeMeasurement(:)]);
            elementsMeasuredIndex = elementsMeasuredIndex + 1;
        end
    end
end
% run through all element names that are to be measured and check if they
% are present in the model.
stateNames = SBstates(model);
variableNames = SBvariables(model);
reactionNames = SBreactions(model);
allNames = {stateNames{:},variableNames{:},reactionNames{:}};
for k1 = 1:length(elementsMeasured),
    foundName = 0;
    for k2 = 1:length(allNames),
        if strcmp(elementsMeasured(k1).name,allNames{k2}),
            foundName = 1;
            break;
        end
    end
    if foundName == 0,
        error('Measurement ''%s'' is not an element of the model.',elementsMeasured(k1).name);
    end
end
% check if any measurements present
if length(elementsMeasured) == 0,
    disp('Nothing selected to measure.');
    data = [];
    return
end
% order the masterTimeVector in ascending order
masterTimeVector = sort(masterTimeVector);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY INFORMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
text = sprintf('Elements selected for measurement:\n');
for k = 1:length(elementsMeasured),
    text = sprintf('%s%s, ',text,elementsMeasured(k).name);
    if mod(k,6) == 0,
        text = sprintf('%s\n',text);
    end    
end
text = text(1:length(text)-2);
disp(text);
notesText = sprintf('%s\n%s',notesText,text);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS THE PERTURBATION INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Processing the desired perturbations to be applied to the model.
% Perturbations are only affecting parameters in the model and are
% implemented using the piecewiseSB function as variables. The
% corresponding parameters need not to be deleted since variables are
% evaluated after the parameters.
experimentmodelStruct = SBstruct(model);
[parameterNames,parameterValuesNominal] = SBparameters(model);
allPerturbationParameters = {};
% check that correct number of elements
if mod(length(perturbations),3) ~= 0,
    error('Incorrect number of elements in the perturbations input argument.');
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY INFORMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
text = 'Perturbations applied to the system:';
disp(text);
notesText = sprintf('%s\n%s',notesText,text);
perturbationParameterNames = {};
% cycle through all the perturbations and update the experiment model structure
for k1 = 1:length(perturbations)/3,
    perturbationName = strtrim(perturbations{3*k1-2});
    perturbationParameterNames{k1} = perturbationName;
    perturbationType = strtrim(perturbations{3*k1-1});
    perturbationValues = perturbations{3*k1};
    % check if perturbationName is the name of a parameter in the model
    foundName = 0;
    for k2 = 1:length(parameterNames),
        if strcmp(perturbationName,parameterNames{k2}),
            foundName = 1;
            parameterIndex = k2;
            break;
        end
    end
    if foundName == 0,
        error('Perturbation parameter ''%s'' not present in model.',perturbationName);
    end
    % check which perturbation type has been chosen, check the format of
    % values and update the model structure for the experiment model by
    % adding a variable for each perturbation
    if strcmp(lower(perturbationType), 'step'),
        % if only 2 values present then add a 1 to set relative
        % perturbations as default perturbations
        if length(perturbationValues) == 2,
            perturbationValues = [perturbationValues(:)' 1];
        end
        % check the length of the values vector
        if length(perturbationValues) ~= 3,
            error('Perturbation values wrongly defined for the perturbation of ''%s''.',perturbationName);
        end
        timeInstant = perturbationValues(1);
        sizePerturbation = perturbationValues(2);
        absRel = perturbationValues(3);
        % process absolute or relative perturbation
        if absRel ~= 0,
            % absRel not zero => do a relative perturbation and assume
            % perturbation size given in percent
            newParameterValue = parameterValuesNominal(parameterIndex)*(1+sizePerturbation/100);
        else
            % absRel zero => do an absoute perturbation
            newParameterValue = parameterValuesNominal(parameterIndex) + sizePerturbation;
        end
        % construct piecewise expression
        variableFormula = sprintf('piecewiseSB(%g,lt(time,%g),%g)',parameterValuesNominal(parameterIndex),timeInstant,newParameterValue);
    elseif strcmp(lower(perturbationType), 'pulse'),
        % if only 3 values present then add a 1 to set relative
        % perturbations as default perturbations
        if length(perturbationValues) == 3,
            perturbationValues = [perturbationValues(:)' 1];
        end
        % check the length of the values vector
        if length(perturbationValues) ~= 4,
            error('Perturbation values wrongly defined for the perturbation of ''%s''.',perturbationName);
        end
        timeInstantChange = perturbationValues(1);
        timeInstantChangeBack = perturbationValues(2);
        % check that change back instant after change instant
        if (timeInstantChange >= timeInstantChangeBack),
            error('The time instants for the parameter change in the perturbation of ''%s'' are wrongly defined.',perturbationName);
        end
        sizePerturbation = perturbationValues(3);
        absRel = perturbationValues(4);
        % process absolute or relative perturbation
        if absRel ~= 0,
            % absRel not zero => do a relative perturbation and assume
            % perturbation size given in percent
            newParameterValue = parameterValuesNominal(parameterIndex)*(1+sizePerturbation/100);
        else
            % absRel zero => do an absoute perturbation
            newParameterValue = parameterValuesNominal(parameterIndex) + sizePerturbation;
        end
        % construct piecewise expression
        variableFormula = sprintf('piecewiseSB(%g,orSB(lt(time,%g),gt(time,%g)),%g)',parameterValuesNominal(parameterIndex),timeInstantChange,timeInstantChangeBack,newParameterValue);
    else
        error('Unknown perturbation type for parameter ''%s''.',perturbationName);
    end
    % now add a new variable of name 'perturbationName' and give it the
    % formula defined by 'variableFormula'.
    % the new variable needs to be added at the top of the variables in the
    % original model. therfore copy the other ones on index higher. 
    variableNumber = length(experimentmodelStruct.variables);
    for k2 = variableNumber:-1:1,
        experimentmodelStruct.variables(k2+1).name = experimentmodelStruct.variables(k2).name;
        experimentmodelStruct.variables(k2+1).formula = experimentmodelStruct.variables(k2).formula;
        experimentmodelStruct.variables(k2+1).type = experimentmodelStruct.variables(k2).type;
        experimentmodelStruct.variables(k2+1).compartment = experimentmodelStruct.variables(k2).compartment;
        experimentmodelStruct.variables(k2+1).unittype = experimentmodelStruct.variables(k2).unittype;
        experimentmodelStruct.variables(k2+1).notes = experimentmodelStruct.variables(k2).notes;
    end
    % add new variable at first place
    experimentmodelStruct.variables(1).name = perturbationName;
    experimentmodelStruct.variables(1).formula = variableFormula;
    experimentmodelStruct.variables(1).type = '';
    experimentmodelStruct.variables(1).compartment = '';
    experimentmodelStruct.variables(1).unittype = '';
    experimentmodelStruct.variables(1).notes = 'Experiment settings set by SBexperiment';
    
    % remember all parameters that are perturbed
    allPerturbationParameters{k1} = perturbationName;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % DISPLAY INFORMATION
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(lower(perturbationType), 'step'),
        text = sprintf('Step perturbation of parameter ''%s''. Nominal value: %g, perturbed value: %g, time of change: %g.',perturbationName,parameterValuesNominal(parameterIndex),newParameterValue,perturbationValues(1));
        notesText = sprintf('%s\n%s',notesText,text);
        disp(text);
    elseif strcmp(lower(perturbationType), 'pulse'),
        text = sprintf('Pulse perturbation of parameter ''%s''. Nominal value: %g, perturbed value: %g, time intervall for changed value: [%g %g].',perturbationName,parameterValuesNominal(parameterIndex),newParameterValue,perturbationValues(1),perturbationValues(2));
        notesText = sprintf('%s\n%s',notesText,text);
        disp(text);
    end
end
% check that no parameter name has been used more than once as perturbation
for k1 = 1:length(allPerturbationParameters),
    for k2 = 1:length(allPerturbationParameters),
        if k1 ~= k2 && strcmp(allPerturbationParameters{k1},allPerturbationParameters{k2}),
            error('Parameter ''%s'' used for more than one perturbation.',allPerturbationParameters{k1});
        end
    end
end
% construct the final experiment model
experimentmodel = SBmodel(experimentmodelStruct);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO THE SIMULATION TO OBTAIN ALL NEEDED DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% simulate the experimental system
output = SBsimulate(experimentmodel,integratorMethod,masterTimeVector,SBinitialconditions(experimentmodel),integratorOptions);
% put all simulation output elements together rowwise
datanames = {output.states{:},output.variables{:},output.reactions{:}};
datavalues = [output.statevalues, output.variablevalues, output.reactionvalues];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS SIMULATION / MEASUREMENT DATA 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% keeping only the measured elements in the right order
elementsMeasureddata = [];
for k = 1:length(elementsMeasured),
    % find the index of the current measurement in the simulation output
    indexMeasurement = 0;
    for k2 = 1:length(datanames),
        if strcmp(elementsMeasured(k).name,datanames{k2}),
            indexMeasurement = k2;
            break;
        end
    end
    if indexMeasurement == 0,
        error('Could not find ''%s'' in simulation output.',elementsMeasured(k).name);
    end
    % get the measurement data 
    elementsMeasureddata(:,k) = datavalues(:,indexMeasurement);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS SIMULATION / STIMULI DATA 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% keeping only the stimulated elements in the right order
elementsStimuliData = [];
for k = 1:length(perturbationParameterNames),
    % find the index of the current measurement in the simulation output
    indexStimulus = 0;
    for k2 = 1:length(datanames),
        if strcmp(datanames{k2},perturbationParameterNames{k}),
            indexStimulus = k2;
            break;
        end
    end
    % get the stimuli data 
    elementsStimuliData(:,k) = datavalues(:,indexStimulus);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS THE NOISE INFORMATION IF PRESENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elementsMeasuredataNoise = [];
if ~isempty(noise),
    % process the noise information
    % first check if noise information correct
    if length(noise) == 2,
        noiseoffset = noise{1};
        noisevariance = noise{2};
        noiseabsRel = 1;
    elseif length(noise) == 3,
        noiseoffset = noise{1};
        noisevariance = noise{2};
        noiseabsRel = noise{3};
    else
        error('Noise input arguments are not correctly given.');
    end
    % check if information given as scalars or if the lengths of the
    % vectors are correct 
    if length(noiseoffset) == 1,
        % offset same for all measurements
        noiseoffset = noiseoffset*ones(1,length(elementsMeasured));
    else
        if length(noiseoffset) ~= length(elementsMeasured),
            error('Size of noise offset vector does not match number of measurements.');
        end
    end
    if length(noisevariance) == 1,
        % variance same for all measurements
        noisevariance = noisevariance*ones(1,length(elementsMeasured));
    else
        if length(noisevariance) ~= length(elementsMeasured),
            error('Size of noise variance vector does not match number of measurements.');
        end
    end
    % now noisevariance and noiseoffset are vectors of length 'number of
    % measurements'
    % adjust the variances and offsets in the case of relative noise data
    if noiseabsRel > 0,
        noiseoffset = noiseoffset(:)'/100.*mean(abs(elementsMeasureddata));
        noisevariance = noisevariance(:)'/100.*mean(abs(elementsMeasureddata));
    end
    % now add the noise to the measurements
    for k = 1:length(elementsMeasured),
        elementsMeasureddataNoise(:,k) = elementsMeasureddata(:,k) + noiseoffset(k) + sqrt(noisevariance(k))*randn(length(masterTimeVector),1);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS THE SIMULATION DATA TO CONSTRUCT SBdata OBJECT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% create empty SBdata object
data = SBdata();
datastructure = SBstruct(data);
datastructurenoise = datastructure;
datastructure.name = strcat('Experiment_on_model_',experimentmodelStruct.name);
datastructurenoise.name = strcat('Experiment_with_noise_on_model_',experimentmodelStruct.name);
datastructure.notes = notesText;
datastructurenoise.notes = notesText;
datastructure.time = masterTimeVector;
datastructurenoise.time = masterTimeVector;
% add measurements
for k = 1:length(elementsMeasured),
    datastructure.measurements(k).name = elementsMeasured(k).name;
    datastructurenoise.measurements(k).name = elementsMeasured(k).name;
    % stimulus flag
    datastructure.measurements(k).stimulusflag = 0;
    datastructurenoise.measurements(k).stimulusflag = 0;
    % fill in noise information
    datastructure.measurements(k).noiseoffset = 0;
    datastructure.measurements(k).noisevariance = 0;
    datastructurenoise.measurements(k).noiseoffset = noiseoffset(k);
    datastructurenoise.measurements(k).noisevariance = noisevariance(k);
    % get time instants to collect for current measurement
    timeIndices = find(ismember(masterTimeVector,elementsMeasured(k).time));
    % get the measurement data 
    datavector = NaN*ones(length(masterTimeVector),1);
    datavector(timeIndices) = elementsMeasureddata(timeIndices,k);
    datastructure.measurements(k).data = datavector(:);
    datavectornoise = NaN*ones(length(masterTimeVector),1);
    datavectornoise(timeIndices) = elementsMeasureddataNoise(timeIndices,k);
    datastructurenoise.measurements(k).data = datavectornoise(:);
end
% add stimuli
for k = 1:length(perturbationParameterNames),
    datastructure.measurements(k+length(elementsMeasured)).name = perturbationParameterNames{k};
    datastructurenoise.measurements(k+length(elementsMeasured)).name = perturbationParameterNames{k};
    % stimulus flag
    datastructure.measurements(k+length(elementsMeasured)).stimulusflag = 1;
    datastructurenoise.measurements(k+length(elementsMeasured)).stimulusflag = 1;
    % fill in noise information
    datastructure.measurements(k+length(elementsMeasured)).noiseoffset = 0;
    datastructure.measurements(k+length(elementsMeasured)).noisevariance = 0;
    datastructurenoise.measurements(k+length(elementsMeasured)).noiseoffset = 0;
    datastructurenoise.measurements(k+length(elementsMeasured)).noisevariance = 0;
    % get the measurement data 
    datastructure.measurements(k+length(elementsMeasured)).data = elementsStimuliData(:,k);
    datastructurenoise.measurements(k+length(elementsMeasured)).data = elementsStimuliData(:,k);
end
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAKE THE DATA OBJECTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data = SBdata(datastructure);
datanoise = SBdata(datastructurenoise);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS CASE OF NO OUTPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout == 0,
    SBvisualizedata(datanoise);
end



return
