function [] = SBexportSBML(sbm, varargin)
% SBexportSBML 
% exports an SBmodel to an SBML Level 2 model using the OutputSBML function
% from the SBML toolbox.
%
% IMPORTANT:
% ==========
% To export a model, the user is required to provide correct information in 
% the "type,compartment,unittype" fields of the internal SBmodel data structure.
%
% This information can be added by directly editing the models structure,
% or by providing the information in the model text representation.
%
% 1) Editing the structure
% ------------------------
% The model structure is obtained by typing "modelstruct = SBstruct(sbmodel)"
% The 'states', 'parameters', and 'variables' fields of the structure have the 
% following subfields that need to be set to correct values.
%   subfields       needed content          comment
%   'type':      	'isSpecie'              in case the state/parameter/variable in
%                                           the SBmodel is an SBML specie
%                   'isParameter'           in case the state/parameter/variable in
%                                           the SBmodel is an SBML parameter
%                   'isCompartment'         in case the state/parameter/variable in
%                                           the SBmodel is an SBML compartment
%
%   'compartment':  specie compartment      if type='isSpecie' it defines the name of the
%                                           compartment the specie is located in.
%                   outside compartment     if type='isCompartment' it defined the name 
%                                           of the compartment outside
%                   ''                      if type='isParameter' or type='isCompartment' 
%                                           and there is no outside compartment or 
%                                           type='isSpecie' and there is no
%                                           compartment defined.
%
%   'unittype':     'amount'                if type='isSpecie' and the species units 
%                                           are in 'amount'
%                   'concentration'         if type='isSpecie' and the species units 
%                                           are in 'concentration'
%                   ''                      if type is not 'isSpecie'                   
%
% 2) Providing the information in the model text description, by SBedit,
%    SBeditBC or directly by editing text files
% ----------------------------------------------------------------------
% For each state, parameter, and variable the additional information has
% the following syntax:
%
% [type:compartment:unittype]
%
% where 'type', 'compartment', and 'unittype' is defined as above.
% This argument needs to be placed in each definition of a state,
% parameter, and variable just before the optional comment.
%
% Example:
% --------
% An example model that has been correctly prepared for SBML export 
% is defined in the file SBTOOLBOX/examples/exampleSBMLexport.txt
%
% USAGE:
% ======
% [] = SBexportSBML(model)
%
% model: SBmodel to export to SBML

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

skipSBMLpreset = 0;
if nargin == 2,
    skipSBMLpreset = varargin{1};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK IF MODEL CONTAINS MATLAB FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(sbm.functionsMATLAB),
    error('Model contains MATLAB functions. Such functions are not supported in SBML.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPROCESS and CHECK SBmodel (fill type, compartment, unittype fields if possible)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~skipSBMLpreset,
    sbm = SBmakeSBMLpresets(sbm);
end
% test wether all name fields are filled and are unique
namesOK = checkNamesInSBmodel(sbm, 1);
if ~namesOK,
    errorMessage{1} = 'SBML export aborted because of errors concerning the names used in the given SBmodel. For further information use "checkNamesInSBmodel" function!';
    messageOutput(errorMessage, 1);
    return    
end
% if test shows that model isn't SBML compliant abort export
[sbmlCompliant, onlyCompartmentError] = testSBMLsettings(sbm, 1);
if ~sbmlCompliant,
    if onlyCompartmentError,
        sbm = addDefaultCompartment(sbm);
    else
        errorMessage{1} = 'SBML export aborted because of errors concerning the SBML compliance. For further information use "testSBMLsettings" function!';
        messageOutput(errorMessage, 1);
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DETERMINE STOICHIOMETRIX MATRIX and related stuff (to determine if
% boundary species or species)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[StoichiometricMatrix,SBMLspeciesNames,SBMLreactionNames,SBMLreversibleFlag] = SBstoichiometry(sbm,1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE AN EMPTY SBML MATLAB STRUCTURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create empty SBMLstructure
% functionDefinition substructure
functionDefinitionStruct = struct('typecode', {}, 'notes', {}, 'annotation', {}, 'name', {}, 'id', {}, 'math', {});
% unitDefinition substructure
unitDefinitionStruct = struct('typecode',{},'notes',{},'annotation',{}, 'name', {}, 'id', {}, 'unit', {});
% compartment substructure
compartmentStruct = struct('typecode', {}, 'notes', {}, 'annotation', {}, 'name', {}, 'id', {}, 'spatialDimensions', {}, 'size', {}, 'units', {}, 'outside', {}, 'constant', {}, 'isSetSize', {}, 'isSetVolume', {});
% species substructure
speciesStruct = struct('typecode', {}, 'notes', {}, 'annotation', {}, 'name', {},'id', {}, 'compartment', {}, 'initialAmount', {}, 'initialConcentration', {}, 'substanceUnits', {}, 'spatialSizeUnits', {}, 'hasOnlySubstanceUnits', {}, 'boundaryCondition', {}, 'charge', {}, 'constant', {}, 'isSetInitialAmount', {}, 'isSetInitialConcentration', {}, 'isSetCharge', {});
% parameter substructure
parameterStruct = struct('typecode', {}, 'notes', {}, 'annotation', {}, 'name', {}, 'id', {}, 'value', {}, 'units', {}, 'constant', {}, 'isSetValue', {});
% rule substructure
ruleStruct = struct('typecode', {}, 'notes', {}, 'annotation', {}, 'formula', {}, 'variable', {}, 'species', {}, 'compartment', {}, 'name', {}, 'units', {});
% reaction substructure
reactantStruct = struct('typecode', {}, 'notes', {}, 'annotation', {}, 'species', {}, 'stoichiometry', {}, 'denominator', {}, 'stoichiometryMath', {});
productStruct = struct('typecode', {}, 'notes', {}, 'annotation', {}, 'species', {}, 'stoichiometry', {}, 'denominator', {}, 'stoichiometryMath', {});
modifierStruct = struct('typecode', {}, 'notes', {}, 'annotation', {}, 'species', {});
kineticLawParameterStruct = struct('typecode', {}, 'notes', {}, 'annotation', {}, 'name', {}, 'id', {}, 'value', {}, 'units', {}, 'constant', {}, 'isSetValue', {});
kineticLawStruct = struct('typecode', {}, 'notes', {}, 'annotation', {}, 'formula', {}, 'math', {}, 'parameter', kineticLawParameterStruct, 'timeUnits', {}, 'substanceUnits', {});
reactionStruct = struct('typecode', {}, 'notes', {}, 'annotation', {}, 'name', {}, 'id', {}, 'reactant', reactantStruct, 'product', productStruct, 'modifier', modifierStruct, 'kineticLaw', kineticLawStruct, 'reversible', {}, 'fast', {}, 'IsSetFast', {});
eventStruct = struct('typecode', {}, 'notes', {}, 'annotation', {}, 'name', {}, 'id', {}, 'trigger', {}, 'delay', {}, 'timeUnits', {}, 'eventAssignment', {});
% Create SBstructure
SBMLstructure = struct('typecode', 'SBML_MODEL',...
                       'notes', '',...
                       'annotation', '', ...
                       'SBML_level', 2, ...
                       'SBML_version', 1, ...
                       'name', '', ...
                       'id', '', ...
                       'functionDefinition', functionDefinitionStruct, ...
                       'unitDefinition', unitDefinitionStruct, ...
                       'compartment', compartmentStruct, ...
                       'species', speciesStruct, ...
                       'parameter', parameterStruct, ...
                       'rule', ruleStruct, ...
                       'reaction', reactionStruct, ...
                       'event', eventStruct, ...
                       'time_symbol', 'time');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOME VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                       
ruleIndex = 1;
parameterIndex = 1;
compartmentIndex = 1;
specieIndex = 1;
                       
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COPY NAME AND NOTES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                       
SBMLstructure.name = sbm.name;
SBMLstructure.notes = convert2SBMLNotes(sbm.notes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS COMPARTMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                       
% compartments can be defined as parameters, variables or states. 
% If the field 'type' is set to 'isCompartment' it indicates a compartment. 
% check parameters for compartments
for k = 1:length(sbm.parameters),
    if strcmp(sbm.parameters(k).type,'isCompartment'),
        % current parameter determines compartment size
        % always use spatial dimensions 3
        SBMLstructure.compartment(compartmentIndex).typecode = 'SBML_COMPARTMENT';
        SBMLstructure.compartment(compartmentIndex).notes = convert2SBMLNotes(sbm.parameters(k).notes);
        SBMLstructure.compartment(compartmentIndex).annotation = '';
        SBMLstructure.compartment(compartmentIndex).name = sbm.parameters(k).name;
        SBMLstructure.compartment(compartmentIndex).id = sbm.parameters(k).name;
        SBMLstructure.compartment(compartmentIndex).spatialDimensions = 3;
        SBMLstructure.compartment(compartmentIndex).size = sbm.parameters(k).value;
        SBMLstructure.compartment(compartmentIndex).units = '';
        SBMLstructure.compartment(compartmentIndex).outside = sbm.parameters(k).compartment;
        SBMLstructure.compartment(compartmentIndex).constant = 1;
        SBMLstructure.compartment(compartmentIndex).isSetSize = 1;
        SBMLstructure.compartment(compartmentIndex).isSetVolume = 1;
        compartmentIndex = compartmentIndex + 1;
    end
end
% check variables for compartments
for k = 1:length(sbm.variables),
    if strcmp(sbm.variables(k).type,'isCompartment'),
        % current variable determines compartment size
        % set size to 1 (not important, since not used)
        SBMLstructure.compartment(compartmentIndex).typecode = 'SBML_COMPARTMENT';
        SBMLstructure.compartment(compartmentIndex).notes = convert2SBMLNotes(sbm.variables(k).notes);
        SBMLstructure.compartment(compartmentIndex).annotation = '';
        SBMLstructure.compartment(compartmentIndex).name = sbm.variables(k).name;
        SBMLstructure.compartment(compartmentIndex).id = sbm.variables(k).name;
        SBMLstructure.compartment(compartmentIndex).spatialDimensions = 3;
        SBMLstructure.compartment(compartmentIndex).size = 1;
        SBMLstructure.compartment(compartmentIndex).units = '';
        SBMLstructure.compartment(compartmentIndex).outside = sbm.variables(k).compartment;
        SBMLstructure.compartment(compartmentIndex).constant = 0;
        SBMLstructure.compartment(compartmentIndex).isSetSize = 1;
        SBMLstructure.compartment(compartmentIndex).isSetVolume = 1;
        compartmentIndex = compartmentIndex + 1;
        % Compartment defined by a rule. Include this rule in the structure
        SBMLstructure.rule(ruleIndex).typecode = 'SBML_ASSIGNMENT_RULE';
        SBMLstructure.rule(ruleIndex).notes = convert2SBMLNotes(sbm.variables(k).notes);
        SBMLstructure.rule(ruleIndex).annotation = '';
        SBMLstructure.rule(ruleIndex).formula = replaceMATLABexpressions(sbm.variables(k).formula);
        SBMLstructure.rule(ruleIndex).variable = sbm.variables(k).name;
        SBMLstructure.rule(ruleIndex).species = '';
        SBMLstructure.rule(ruleIndex).compartment = '';
        SBMLstructure.rule(ruleIndex).name = '';
        SBMLstructure.rule(ruleIndex).units = '';
        ruleIndex = ruleIndex + 1;
    end
end
% check states for compartments
for k = 1:length(sbm.states),
    if strcmp(sbm.states(k).type,'isCompartment'),
        % current state determines compartment size
        % set size to initial condition of state
        SBMLstructure.compartment(compartmentIndex).typecode = 'SBML_COMPARTMENT';
        SBMLstructure.compartment(compartmentIndex).notes = convert2SBMLNotes(sbm.states(k).notes);
        SBMLstructure.compartment(compartmentIndex).annotation = '';
        SBMLstructure.compartment(compartmentIndex).name = sbm.states(k).name;
        SBMLstructure.compartment(compartmentIndex).id = sbm.states(k).name;
        SBMLstructure.compartment(compartmentIndex).spatialDimensions = 3;
        SBMLstructure.compartment(compartmentIndex).size = sbm.states(k).initialCondition;
        SBMLstructure.compartment(compartmentIndex).units = '';
        SBMLstructure.compartment(compartmentIndex).outside = sbm.states(k).compartment;
        SBMLstructure.compartment(compartmentIndex).constant = 0;
        SBMLstructure.compartment(compartmentIndex).isSetSize = 1;
        SBMLstructure.compartment(compartmentIndex).isSetVolume = 1;
        compartmentIndex = compartmentIndex + 1;
        % Compartment defined by a rule. Include this rule in the structure
        SBMLstructure.rule(ruleIndex).typecode = 'SBML_RATE_RULE';
        SBMLstructure.rule(ruleIndex).notes = convert2SBMLNotes(sbm.states(k).notes);
        SBMLstructure.rule(ruleIndex).annotation = '';
        SBMLstructure.rule(ruleIndex).formula = replaceMATLABexpressions(sbm.states(k).ODE);
        SBMLstructure.rule(ruleIndex).variable = sbm.states(k).name;
        SBMLstructure.rule(ruleIndex).species = '';
        SBMLstructure.rule(ruleIndex).compartment = '';
        SBMLstructure.rule(ruleIndex).name = '';
        SBMLstructure.rule(ruleIndex).units = '';
        ruleIndex = ruleIndex + 1;        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% parameters can be defined as parameters, variables or states. If the
% field 'type' is set to 'isParameter' it indicates an SBML parameter.  
% we only will have global parameters in the SBML model.
% check parameters for parameters
for k = 1:length(sbm.parameters),
    if strcmp(sbm.parameters(k).type,'isParameter'),
        % current parameter determines compartment size
        SBMLstructure.parameter(parameterIndex).typecode = 'SBML_PARAMETER';
        SBMLstructure.parameter(parameterIndex).notes = convert2SBMLNotes(sbm.parameters(k).notes);
        SBMLstructure.parameter(parameterIndex).annotation = '';
        SBMLstructure.parameter(parameterIndex).name = sbm.parameters(k).name;
        SBMLstructure.parameter(parameterIndex).id = sbm.parameters(k).name;
        SBMLstructure.parameter(parameterIndex).value = sbm.parameters(k).value;
        SBMLstructure.parameter(parameterIndex).units = '';
        SBMLstructure.parameter(parameterIndex).constant = 1;
        SBMLstructure.parameter(parameterIndex).isSetValue = 1;
        parameterIndex = parameterIndex + 1;
    end
end
% check variables for parameters
for k = 1:length(sbm.variables),
    if strcmp(sbm.variables(k).type,'isParameter'),
        % current variable determines compartment size
        % value set to zero (not needed)
        SBMLstructure.parameter(parameterIndex).typecode = 'SBML_PARAMETER';
        SBMLstructure.parameter(parameterIndex).notes = convert2SBMLNotes(sbm.variables(k).notes);
        SBMLstructure.parameter(parameterIndex).annotation = '';
        SBMLstructure.parameter(parameterIndex).name = sbm.variables(k).name;
        SBMLstructure.parameter(parameterIndex).id = sbm.variables(k).name;
        SBMLstructure.parameter(parameterIndex).value = 0;  % no value needed (assignment rule)
        SBMLstructure.parameter(parameterIndex).units = '';
        SBMLstructure.parameter(parameterIndex).constant = 0;
        SBMLstructure.parameter(parameterIndex).isSetValue = 1;
        parameterIndex = parameterIndex + 1;        
        % Parameter defined by a rule. Include this rule in the structure
        SBMLstructure.rule(ruleIndex).typecode = 'SBML_ASSIGNMENT_RULE';
        SBMLstructure.rule(ruleIndex).notes = convert2SBMLNotes(sbm.variables(k).notes);
        SBMLstructure.rule(ruleIndex).annotation = '';
        SBMLstructure.rule(ruleIndex).formula = replaceMATLABexpressions(sbm.variables(k).formula);
        SBMLstructure.rule(ruleIndex).variable = sbm.variables(k).name;
        SBMLstructure.rule(ruleIndex).species = '';
        SBMLstructure.rule(ruleIndex).compartment = '';
        SBMLstructure.rule(ruleIndex).name = '';
        SBMLstructure.rule(ruleIndex).units = '';
        ruleIndex = ruleIndex + 1;
    end
end
% check states for parameters
for k = 1:length(sbm.states),
    if strcmp(sbm.states(k).type,'isParameter'),
        % current state determines compartment size
        SBMLstructure.parameter(parameterIndex).typecode = 'SBML_PARAMETER';
        SBMLstructure.parameter(parameterIndex).notes = convert2SBMLNotes(sbm.states(k).notes);
        SBMLstructure.parameter(parameterIndex).annotation = '';
        SBMLstructure.parameter(parameterIndex).name = sbm.states(k).name;
        SBMLstructure.parameter(parameterIndex).id = sbm.states(k).name;
        SBMLstructure.parameter(parameterIndex).value = sbm.states(k).initialCondition;
        SBMLstructure.parameter(parameterIndex).units = '';
        SBMLstructure.parameter(parameterIndex).constant = 0;
        SBMLstructure.parameter(parameterIndex).isSetValue = 1;
        parameterIndex = parameterIndex + 1;        
        % Parameter defined by a rule. Include this rule in the structure
        SBMLstructure.rule(ruleIndex).typecode = 'SBML_RATE_RULE';
        SBMLstructure.rule(ruleIndex).notes = convert2SBMLNotes(sbm.states(k).notes);
        SBMLstructure.rule(ruleIndex).annotation = '';
        SBMLstructure.rule(ruleIndex).formula = replaceMATLABexpressions(sbm.states(k).ODE);
        SBMLstructure.rule(ruleIndex).variable = sbm.states(k).name;
        SBMLstructure.rule(ruleIndex).species = '';
        SBMLstructure.rule(ruleIndex).compartment = '';
        SBMLstructure.rule(ruleIndex).name = '';
        SBMLstructure.rule(ruleIndex).units = '';
        ruleIndex = ruleIndex + 1;        
    end
end            

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% add functions
for k = 1:length(sbm.functions),
    SBMLstructure.functionDefinition(k).typecode = 'SBML_FUNCTION_DEFINITION';
    SBMLstructure.functionDefinition(k).notes = convert2SBMLNotes(sbm.functions(k).notes);
    SBMLstructure.functionDefinition(k).annotation = '';
    SBMLstructure.functionDefinition(k).name = sbm.functions(k).name;
    SBMLstructure.functionDefinition(k).id = sbm.functions(k).name;
    SBMLstructure.functionDefinition(k).math = sprintf('lambda(%s,%s)',sbm.functions(k).arguments,replaceMATLABexpressions(sbm.functions(k).formula));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS CONSTANT (BOUNDARY) SPECIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% constant species can only be defined as parameters.
% If the field 'type' of a parameter is set to 'isSpecie' we assume that
% this indicates a constant boundary species.
for k = 1:length(sbm.parameters),
    if strcmp(sbm.parameters(k).type,'isSpecie'),
        % Determine if amount or concentration
        if strcmp(sbm.parameters(k).unittype,'amount') 
            isAmount = 1;
            isConcentration = 0;
        elseif strcmp(sbm.parameters(k).unittype,'concentration') 
            isAmount = 0;
            isConcentration = 1;
        end
        % import constant boundary species
        SBMLstructure.species(specieIndex).typecode = 'SBML_SPECIES';
        SBMLstructure.species(specieIndex).notes = convert2SBMLNotes(sbm.parameters(k).notes);
        SBMLstructure.species(specieIndex).annotation = '';
        SBMLstructure.species(specieIndex).name = sbm.parameters(k).name;
        SBMLstructure.species(specieIndex).id = sbm.parameters(k).name;
        SBMLstructure.species(specieIndex).compartment = sbm.parameters(k).compartment;
        SBMLstructure.species(specieIndex).initialAmount = sbm.parameters(k).value;
        SBMLstructure.species(specieIndex).initialConcentration = sbm.parameters(k).value;
        SBMLstructure.species(specieIndex).substanceUnits = '';         % no units!
        SBMLstructure.species(specieIndex).spatialSizeUnits = '';       % no units!
        SBMLstructure.species(specieIndex).hasOnlySubstanceUnits = 0;
        SBMLstructure.species(specieIndex).boundaryCondition = 1;       % boundary 
        SBMLstructure.species(specieIndex).charge = 0;                  
        SBMLstructure.species(specieIndex).constant = 1;                % species is set constant
        SBMLstructure.species(specieIndex).isSetInitialAmount = isAmount;
        SBMLstructure.species(specieIndex).isSetInitialConcentration = isConcentration;
        SBMLstructure.species(specieIndex).isSetCharge = 0;             % charge not supported
        specieIndex = specieIndex + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS BOUNDARY SPECIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% boundary species can be defined variables or states (not as parameters,
% since then they become constant boundary species).
% The difference between species and boundary species is here that species
% are only defined by reactions, while boundary species are defined by
% rules. The variable "SBMLspeciesNames" contains the names of all states
% that are defined only by reactions.
% now check variables (all species in variables are boundary species)
for k = 1:length(sbm.variables),
    if strcmp(sbm.variables(k).type,'isSpecie'),
        % Determine if amount or concentration
        if strcmp(sbm.variables(k).unittype,'amount') 
            isAmount = 1;
            isConcentration = 0;
        elseif strcmp(sbm.variables(k).unittype,'concentration') 
            isAmount = 0;
            isConcentration = 1;
        end
        % import boundary species
        SBMLstructure.species(specieIndex).typecode = 'SBML_SPECIES';
        SBMLstructure.species(specieIndex).notes = convert2SBMLNotes(sbm.variables(k).notes);
        SBMLstructure.species(specieIndex).annotation = '';
        SBMLstructure.species(specieIndex).name = sbm.variables(k).name;
        SBMLstructure.species(specieIndex).id = sbm.variables(k).name;
        SBMLstructure.species(specieIndex).compartment = sbm.variables(k).compartment;
        SBMLstructure.species(specieIndex).initialAmount = 0;  % not set (assignment rule)
        SBMLstructure.species(specieIndex).initialConcentration = 0; % not set (assignment rule)
        SBMLstructure.species(specieIndex).substanceUnits = '';         % no units!
        SBMLstructure.species(specieIndex).spatialSizeUnits = '';       % no units!
        SBMLstructure.species(specieIndex).hasOnlySubstanceUnits = 0;
        SBMLstructure.species(specieIndex).boundaryCondition = 1;       % boundary 
        SBMLstructure.species(specieIndex).charge = 0;                  
        SBMLstructure.species(specieIndex).constant = 0;                % species is not set constant
        SBMLstructure.species(specieIndex).isSetInitialAmount = isAmount;
        SBMLstructure.species(specieIndex).isSetInitialConcentration = isConcentration;
        SBMLstructure.species(specieIndex).isSetCharge = 0;             % charge not supported
        specieIndex = specieIndex + 1;
        % Boundary species defined by a rule. Include this rule in the structure
        SBMLstructure.rule(ruleIndex).typecode = 'SBML_ASSIGNMENT_RULE';
        SBMLstructure.rule(ruleIndex).notes = convert2SBMLNotes(sbm.variables(k).notes);
        SBMLstructure.rule(ruleIndex).annotation = '';
        SBMLstructure.rule(ruleIndex).formula = replaceMATLABexpressions(sbm.variables(k).formula);
        SBMLstructure.rule(ruleIndex).variable = sbm.variables(k).name;
        SBMLstructure.rule(ruleIndex).species = '';
        SBMLstructure.rule(ruleIndex).compartment = '';
        SBMLstructure.rule(ruleIndex).name = '';
        SBMLstructure.rule(ruleIndex).units = '';
        ruleIndex = ruleIndex + 1;
    end
end
% now check states (boundary if isSpecie and not in the SBMLspeciesNames list)
for k = 1:length(sbm.states),
    if strcmp(sbm.states(k).type,'isSpecie') && isempty(strmatch(sbm.states(k).name,SBMLspeciesNames)),
        % Determine if amount or concentration
        if strcmp(sbm.states(k).unittype,'amount') 
            isAmount = 1;
            isConcentration = 0;
        elseif strcmp(sbm.states(k).unittype,'concentration') 
            isAmount = 0;
            isConcentration = 1;
        end
        % import boundary species
        SBMLstructure.species(specieIndex).typecode = 'SBML_SPECIES';
        SBMLstructure.species(specieIndex).notes = convert2SBMLNotes(sbm.states(k).notes);
        SBMLstructure.species(specieIndex).annotation = '';
        SBMLstructure.species(specieIndex).name = sbm.states(k).name;
        SBMLstructure.species(specieIndex).id = sbm.states(k).name;
        SBMLstructure.species(specieIndex).compartment = sbm.states(k).compartment;
        SBMLstructure.species(specieIndex).initialAmount = sbm.states(k).initialCondition;
        SBMLstructure.species(specieIndex).initialConcentration = sbm.states(k).initialCondition;
        SBMLstructure.species(specieIndex).substanceUnits = '';         % no units!
        SBMLstructure.species(specieIndex).spatialSizeUnits = '';       % no units!
        SBMLstructure.species(specieIndex).hasOnlySubstanceUnits = 0;
        SBMLstructure.species(specieIndex).boundaryCondition = 1;       % boundary 
        SBMLstructure.species(specieIndex).charge = 0;                  
        SBMLstructure.species(specieIndex).constant = 0;                % species is not set constant
        SBMLstructure.species(specieIndex).isSetInitialAmount = isAmount;
        SBMLstructure.species(specieIndex).isSetInitialConcentration = isConcentration;
        SBMLstructure.species(specieIndex).isSetCharge = 0;             % charge not supported
        specieIndex = specieIndex + 1;
        % Boundary species defined by a rule. Include this rule in the structure
        SBMLstructure.rule(ruleIndex).typecode = 'SBML_RATE_RULE';
        SBMLstructure.rule(ruleIndex).notes = convert2SBMLNotes(sbm.states(k).notes);
        SBMLstructure.rule(ruleIndex).annotation = '';
        SBMLstructure.rule(ruleIndex).formula = replaceMATLABexpressions(sbm.states(k).ODE);
        SBMLstructure.rule(ruleIndex).variable = sbm.states(k).name;
        SBMLstructure.rule(ruleIndex).species = '';
        SBMLstructure.rule(ruleIndex).compartment = '';
        SBMLstructure.rule(ruleIndex).name = '';
        SBMLstructure.rule(ruleIndex).units = '';
        ruleIndex = ruleIndex + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS SPECIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% we assume that species are only defined as states and determined by
% reactions
% check simple species defined by reactions
% save the species name and the species ODE for use when processing the
% reactions
speciesReactions = [];
for k = 1:length(sbm.states),
    if strcmp(sbm.states(k).type,'isSpecie') && ~isempty(strmatch(sbm.states(k).name,SBMLspeciesNames)),
        % Determine if amount or concentration
        if strcmp(sbm.states(k).unittype,'amount') 
            isAmount = 1;
            isConcentration = 0;
        elseif strcmp(sbm.states(k).unittype,'concentration') 
            isAmount = 0;
            isConcentration = 1;
        end
        % import constant species
        SBMLstructure.species(specieIndex).typecode = 'SBML_SPECIES';
        SBMLstructure.species(specieIndex).notes = convert2SBMLNotes(sbm.states(k).notes);
        SBMLstructure.species(specieIndex).annotation = '';
        SBMLstructure.species(specieIndex).name = sbm.states(k).name;
        SBMLstructure.species(specieIndex).id = sbm.states(k).name;
        SBMLstructure.species(specieIndex).compartment = sbm.states(k).compartment;
        SBMLstructure.species(specieIndex).initialAmount = sbm.states(k).initialCondition;
        SBMLstructure.species(specieIndex).initialConcentration = sbm.states(k).initialCondition;
        SBMLstructure.species(specieIndex).substanceUnits = '';         % no units!
        SBMLstructure.species(specieIndex).spatialSizeUnits = '';       % no units!
        SBMLstructure.species(specieIndex).hasOnlySubstanceUnits = 0;
        SBMLstructure.species(specieIndex).boundaryCondition = 0;       % not boundary 
        SBMLstructure.species(specieIndex).charge = 0;                  
        SBMLstructure.species(specieIndex).constant = 0;                % species is not set constant
        SBMLstructure.species(specieIndex).isSetInitialAmount = isAmount;
        SBMLstructure.species(specieIndex).isSetInitialConcentration = isConcentration;
        SBMLstructure.species(specieIndex).isSetCharge = 0;             % charge not supported
        specieIndex = specieIndex + 1;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS REACTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Before processing the reactions we are checking for moiety conservations
% and reversing them into the stoichiometric matrix (to be able to draw the
% full picture)
[StoichiometricMatrix,SBMLspeciesNames] = reverseMoietyConservationsSB(sbm, 1);
for k1 = 1:length(sbm.reactions),
    reactionName = sbm.reactions(k1).name;
    reactionFormula = sbm.reactions(k1).formula;
    % get needed information about reaction
    [reactant, product, modifier] = getReactionInfo(reactionName,StoichiometricMatrix,SBMLspeciesNames,SBMLreactionNames,SBMLreversibleFlag,sbm);
    % add reaction to the structure
    SBMLstructure.reaction(k1).typecode = 'SBML_REACTION';
    SBMLstructure.reaction(k1).notes = convert2SBMLNotes(sbm.reactions(k1).notes);
    SBMLstructure.reaction(k1).annotation = '';
    SBMLstructure.reaction(k1).name = reactionName;
    SBMLstructure.reaction(k1).id = reactionName;
    
    SBMLstructure.reaction(k1).reactant = reactant;
    SBMLstructure.reaction(k1).product = product;
    SBMLstructure.reaction(k1).modifier = modifier; 
    
    SBMLstructure.reaction(k1).kineticLaw.typecode = 'SBML_KINETIC_LAW';
    SBMLstructure.reaction(k1).kineticLaw.notes = '';
    SBMLstructure.reaction(k1).kineticLaw.annotation = '';
    SBMLstructure.reaction(k1).kineticLaw.formula = replaceMATLABexpressions(reactionFormula);
    SBMLstructure.reaction(k1).kineticLaw.math = replaceMATLABexpressions(reactionFormula);
    SBMLstructure.reaction(k1).kineticLaw.parameter = kineticLawParameterStruct;    % no local parameters
    
    SBMLstructure.reaction(k1).kineticLaw.timeUnits = '';
    SBMLstructure.reaction(k1).kineticLaw.substanceUnits = '';
    
    SBMLstructure.reaction(k1).reversible = sbm.reactions(k1).reversible;
    
    SBMLstructure.reaction(k1).fast = 0;
    SBMLstructure.reaction(k1).IsSetFast = 0;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LEAVE UNITDEFINITIONS EMPTY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process Events
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                       
for k1=1:length(sbm.events),
    SBMLstructure.event(k1).typecode = 'SBML_EVENT';
    SBMLstructure.event(k1).notes = convert2SBMLNotes(sbm.events(k1).notes);
    SBMLstructure.event(k1).annotation = '';
    SBMLstructure.event(k1).name = sbm.events(k1).name;
    SBMLstructure.event(k1).id = sbm.events(k1).name;
    SBMLstructure.event(k1).trigger = replaceMATLABexpressions(sbm.events(k1).trigger);
    SBMLstructure.event(k1).delay = '0'; % always zero
    SBMLstructure.event(k1).timeUnits = '';
    for k2 = 1:length(sbm.events(k1).assignment),
        SBMLstructure.event(k1).eventAssignment(k2).typecode = 'SBML_EVENT_ASSIGNMENT';
        SBMLstructure.event(k1).eventAssignment(k2).notes = '';
        SBMLstructure.event(k1).eventAssignment(k2).annotation = '';
        SBMLstructure.event(k1).eventAssignment(k2).variable = sbm.events(k1).assignment(k2).variable;
        SBMLstructure.event(k1).eventAssignment(k2).math =  replaceMATLABexpressions(sbm.events(k1).assignment(k2).formula);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITE SBML FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                       
result = SBMLstructure;
try 
    OutputSBML(SBMLstructure);
catch
    if ~strcmp(lasterr, 'Cannot copy filename'),
        warning(sprintf('Error occurred: %s\n',lasterr));
    end
end
return



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DETERMINE REACTION INFORMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                       
% we look through all ODEs where the gioven reaction is involved in.
% we determine reversibility, reactants, products, and modifiers.
function [reactant, product, modifier] = getReactionInfo(reactionName,StoichiometricMatrix,SBMLspeciesNames,SBMLreactionNames,SBMLreversibleFlag,sbm)
    % initialize the structures
    reactant = struct('typecode', {}, 'notes', {}, 'annotation', {}, 'species', {}, 'stoichiometry', {}, 'denominator', {}, 'stoichiometryMath', {});
    product = struct('typecode', {}, 'notes', {}, 'annotation', {}, 'species', {}, 'stoichiometry', {}, 'denominator', {}, 'stoichiometryMath', {});
    modifier = struct('typecode', {}, 'notes', {}, 'annotation', {}, 'species', {});
    % get index of current reaction
    reactionIndex = strmatch(reactionName,SBMLreactionNames,'exact');
    % get corresponding stoichiometries
    [stoichRows, stoichCols] = size(StoichiometricMatrix);
    if (stoichCols >= reactionIndex),
        Nreaction = StoichiometricMatrix(:,reactionIndex);
    else
        % for the current reaction there's no stoichiometry
        % or products, reactants or modifierer found
        return
    end
    % get reactant indices (Nreaction < 0)
    reactantIndices = find(Nreaction < 0);
    % get product indices (Nreaction > 0)
    productIndices = find(Nreaction > 0);
    % add reactants
    for k = 1:length(reactantIndices),
        reactant(k).typecode = 'SBML_SPECIES_REFERENCE';
        reactant(k).notes = '';
        reactant(k).annotation = '';
        reactant(k).species = SBMLspeciesNames{reactantIndices(k)};
        reactant(k).stoichiometry = abs(Nreaction(reactantIndices(k)));
        reactant(k).denominator = 1;
        reactant(k).stoichiometryMath = '';
    end
    % add products
    for k = 1:length(productIndices),
        product(k).typecode = 'SBML_SPECIES_REFERENCE';
        product(k).notes = '';
        product(k).annotation = '';
        product(k).species = SBMLspeciesNames{productIndices(k)};
        product(k).stoichiometry = Nreaction(productIndices(k));
        product(k).denominator = 1;
        product(k).stoichiometryMath = '';
    end
    % search for modifiers among the non reactants and non products and add them
    specieListPossible = SBMLspeciesNames(find(Nreaction == 0));
    modifierList = {};
    modifierList = search4Modifiers(modifierList,sbm.reactions(reactionIndex).formula,specieListPossible,{sbm.variables.name},sbm);
    modifierList = unique(modifierList);
    for k = 1:length(modifierList),
        modifier(k).typecode = 'SBML_MODIFIER_SPECIES_REFERENCE';
        modifier(k).notes = '';
        modifier(k).annotation = '';
        modifier(k).species = modifierList{k};
    end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DETERMINE REACTION INFORMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                       
% we go through all reactions and eventually variables to determine the
% modifiers (recursively)
function [modifierList] = search4Modifiers(modifierList,formula,specieListPossible,variableNames,sbm)
    % search in current formula for modifiers
    formula = strcat('#',formula,'#');
    for k=1:length(specieListPossible),
        if ~isempty(regexp(formula,strcat('\W',specieListPossible{k},'\W'))),
            modifierList{end+1} = specieListPossible{k};
        end
    end
    % search in current formula for variable names
    searchVariablesIndices = [];
    for k=1:length(variableNames),
        if ~isempty(regexp(formula,strcat('\W',variableNames{k},'\W'))),
            searchVariablesIndices(end+1) = k;
        end
    end
    % loop through variable names and search for modifiers
    for k=1:length(searchVariablesIndices),
        formula = sbm.variables(searchVariablesIndices(k)).formula;
        modifierList = search4Modifiers(modifierList,formula,specieListPossible,variableNames,sbm);
    end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXCHANGE MATLAB FUNCTION NAMES TO MATHML FUNCTION NAMES IN ALL FORMULAS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MathML functions can have different names than MATLAB function names
% we need to exchange them in the formula strings. The term 'MathML
% expression' relates to what is required by OutputSBML!
function [newFormula] = replaceMATLABexpressions(formula)
% MathML expressions 
MathMLexpressions = {'arccos(','arcsin(','arctan(','ceiling(','ln(','log(10,',...
    'pow(','arccos(','arccosh(','arccot(','arccoth(','arccsc(','arccsch(',...
    'arcsec(','arcsech(','arcsin(','arcsinh(','arctan(','arctanh(',...
    'exponentiale','log('};
% Corresponding MATLAB expressions
MATLABexpressions = {'acos(','asin(','atan(','ceil(','log(','log10(',...
    'power(','acos(','acosh(','acot(','acoth(','acsc(','acsch(',...
    'asec(','asech(','asin(','asinh(','atan(','atanh('...
    'exp(1)','logbase('};
% do the simple replace
for k = 1:length(MathMLexpressions),
    newFormula = strrep(formula,MATLABexpressions{k},MathMLexpressions{k});
end
return


    
    
    