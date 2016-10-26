function [SBstructure,errorMsg] = SBconvertSBML1(SBMLmodel)
% SBconvertSBML1
% Converting a SBML Level 1 MATLAB structure to the object model structure
% used in the toolbox. Works for SBML Level 1, versions 1 and 2
%
% Private method 
%
% [SBstructure,errorMsg] = SBconvertSBML1(SBMLmodel);
%
% Simple error checking is provided. In the case of an error an empty
% structure is returned. 

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
% INITIALIZE ERROR MESSAGES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
errorMsg = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CREATE AN EMPTY SBSTRUCTURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create empty SBstructure
% functions substructure
functionsStruct = struct('name',{},'arguments',{},'formula',{},'notes',{});
% states substructure
statesStruct = struct('name',{},'initialCondition',{},'ODE',{},'type',{},'compartment',{},'unittype',{},'notes',{});
% parameters substructure
parametersStruct = struct('name',{},'value',{},'type',{},'compartment',{},'unittype',{},'notes',{});
% variables substructure
variablesStruct = struct('name',{},'formula',{},'type',{},'compartment',{},'unittype',{},'notes',{});
% reactions substructure
reactionsStruct = struct('name',{},'formula',{},'notes',{},'reversible',{});
% event assignment substructure
eventassignmentStruct = struct('variable',{},'formula',{});
% event substructure
eventStruct = struct('name',{},'trigger',{},'assignment',eventassignmentStruct,'notes',{});
% Create SBstructure
SBstructure = struct('name','','notes','','functions',functionsStruct,'states',statesStruct,'parameters',parametersStruct,'variables',variablesStruct,'reactions',reactionsStruct,'events',eventStruct,'functionsMATLAB','');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK NOT SUPPORTED RULES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Algebraic rules can not be supported and are to avoid. The presence of
% algebraic rules leads to an error. 
if ~isempty(strmatch('SBML_ALGEBRAIC_RULE',{SBMLmodel.rule.typecode})),
    errorMsg = sprintf('%s\nAn algebraic rule is present in the SBML model. These rules are not\nsupported by this toolbox\n',errorMsg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NAME AND NOTES (INCLUDE IN STRUCTURE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(SBMLmodel.name,''),
    % No name provided by the SBML model
    SBstructure.name = 'Unknown name';
else
    % Use the name and the notes supplied in the SBMLmodel
    SBstructure.name = SBMLmodel.name;
end
% Indicate in notes that the model has been created by import from SBML.
SBstructure.notes = SBMLmodel.notes;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOME COUNTER DEFINITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% States and parameters are included in the structure in their order of
% appearance. For variables no counters are needed, since their order 
% is defined by the ordering of the scalar rules and not by their order of 
% appearance (see below)
numberStates = 1;
numberParameters = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ORDERING OF VARIABLES (SCALAR RULES)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Scalar rules lead always to the definition of variables. The order of
% appearance of the rules is important and thus the order of appearance of
% the variables should reflect the ordering of the rules. We determine an
% array in which the the index i determines the priority of the scalar rule, and
% the i-th array element the index of the corresponding scalar rule. Rate
% rules do not have to be considered since they are evaluated at last
% anyway
orderScalarRules = strmatch('SBML_ASSIGNMENT_RULE',{SBMLmodel.rule.typecode});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SPECIES (INCLUDE IN STRUCTURE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Species of the model are included in the structure as follows:
% Non-boundary species:
%   - rate rule: as state (not supported by libSBML yet)
%   - scalar rule: as variable 
%   - no rule: as state
%
% Boundary species:
%   - rate rule: as states (not supported by libSBML yet)
%   - scalar rule: as variable
%   - no rule: as parameter
%
% Initial values for species are assumed to be given in amount units.
% They are converted to concentration units by dividing with the
% compartment size. Specie values set by rules are not converted and taken
% as they are.
%
speciesODElist = [];
speciesODElistIndex = 1;
% Check if there are any species. Otherwise an error message is set
if length(SBMLmodel.species)==0,
%    errorMsg = sprintf('%s\n%s\n',errorMsg,'No species are defined in the SBML model');
    warning('No species are defined in the SBML model - uncertain if model conversion will lead to a working SBmodel!');
end
% Cycle through all the species and determine how they are to be included
% in the model structure
for k1 = 1:length(SBMLmodel.species),
    speciesName = SBMLmodel.species(k1).name;
    initialAmount = SBMLmodel.species(k1).initialAmount;
    speciesCompartmentName = SBMLmodel.species(k1).compartment;
    % get the compartment size for this species
    for k2 = 1:length(SBMLmodel.compartment),
        if strcmp(speciesCompartmentName, SBMLmodel.compartment(k2).name),
            compartmentSize = SBMLmodel.compartment(k2).volume;
            if length(compartmentSize)==0,
                errorMsg = sprintf('%s\nNo compartment size defined for compartment ''%s''\n',errorMsg,speciesCompartmentName);
                compartmentSize = 1; % set to whatever, error will still occurr!
            end
            break;
        end
    end
    % Check the number of rules the species is in and get the index
    % of the last rule the species is used in together with its type
    % ('rate' or 'scalar') and the right hand side expression 'formula'
    [numberRulesSpecies,lastRuleIndex,lastRuleType,lastRuleFormula,errorMsg] = getRules(speciesName,SBMLmodel,errorMsg);
    if numberRulesSpecies > 1,
        % Species defined by more than one rule => error
        errorMsg = sprintf('%s\nSpecies ''%s'' defined by more than one rule\n',errorMsg,speciesName);
    end
    if SBMLmodel.species(k1).boundaryCondition == 0,
        % Process all the species that DO NOT have a boundary condition
        if numberRulesSpecies == 0,
            % add species to the ODE list (for use below in ODE
            % construction)
            speciesODElist(speciesODElistIndex).name = speciesName;
            speciesODElist(speciesODElistIndex).stateIndex = numberStates;
            speciesODElistIndex = speciesODElistIndex + 1;
            % Include species as state
            SBstructure.states(numberStates).name = speciesName;
            % Check if an initial amount is set for a species
            if length(initialAmount)==0,
                errorMsg = sprintf('%s\nNo initial amount defined for species ''%s''\n',errorMsg,speciesName);
            end           
            % Convert initialAmount to concentration and include in SBstructure
            SBstructure.states(numberStates).initialCondition = initialAmount/compartmentSize;
            % Note that this state is a 'species'
            % Notes are not used by the toolbox for other things as
            % documentation
            SBstructure.states(numberStates).notes = '';
            % initialize type, compartment, and unittype fields
            SBstructure.states(numberStates).type = 'isSpecie';
            SBstructure.states(numberStates).compartment = SBMLmodel.species(k1).compartment;
            SBstructure.states(numberStates).unittype = 'concentration';                
            % Initialize the ODE field for this state
            SBstructure.states(numberStates).ODE = '';
            % Increment
            numberStates = numberStates+1;
        elseif numberRulesSpecies == 1,
            % One rule has been detected for this species
            % Before including the species in the model structure we have
            % to check if this species is also used in reactions as product
            % or reactant. In this case an error message will be issued.
            reactionSpeciesPresent = checkReactionsForSpecies(speciesName,SBMLmodel);
            if reactionSpeciesPresent,
                % Species is also present in a reaction. This means the SBML
                % model is inconsistent. => error
                errorMsg = sprintf('%s\nSpecies ''%s'' defined by rule and altered by reactions\n',errorMsg,speciesName);
            end
            % Even in the case of an error we continue here with the
            % processing - The error message will take care of returning
            % an empty model structure at the end
            if strcmp(lastRuleType,'rate'),
                % Species added as a state, since rule of type 'rate'
                SBstructure.states(numberStates).name = speciesName;
                % check if an initial amount is set for a species
                if length(initialAmount)==0,
                    errorMsg = sprintf('%s\nNo initial amount defined for species ''%s''\n',errorMsg,speciesName);
                end           
                % Convert initialAmount to concentration and include in SBstructure
                SBstructure.states(numberStates).initialCondition = initialAmount/compartmentSize;
                % Include the RHS formula as the ODE for this species state
                SBstructure.states(numberStates).ODE = lastRuleFormula;
                % Set a note that type of state is 'species'
                SBstructure.states(numberStates).notes = '';
                % initialize type, compartment, and unittype fields
                SBstructure.states(numberStates).type = 'isParameter';
                SBstructure.states(numberStates).compartment = '';
                SBstructure.states(numberStates).unittype = '';                
                % Increment
                numberStates = numberStates+1;
            else
                % Species added as a variable, since rule of type 'scalar'
                % Determine the index where to include the variable (based
                % on the ordering of the scalar rules)
                indexVariable = find(orderScalarRules==lastRuleIndex);
                SBstructure.variables(indexVariable).name = speciesName;
                SBstructure.variables(indexVariable).formula = lastRuleFormula;
                SBstructure.variables(indexVariable).notes = '';
                % initialize type, compartment, and unittype fields
                SBstructure.variables(indexVariable).type = 'isParameter';
                SBstructure.variables(indexVariable).compartment = '';
                SBstructure.variables(indexVariable).unittype = '';
            end
        end
        % The case where numberRulesSpecies > 1 does not have to be treated 
        % here since already an error message is set above. This error
        % message will lead to an empty model structure at the end. Even 
        % if the error is set the rest of the model will be processed.
        % This has the advantage that all errors can be
        % detected at once and makes the control structure of this script
        % simpler.
    else
        % Process all the species that DO have a boundary condition
        if numberRulesSpecies == 0,
            % Include boundary species as a parameter
            SBstructure.parameters(numberParameters).name = speciesName;
            % Check if an initial amount is set for a species
            if length(initialAmount)==0,
                errorMsg = sprintf('%s\nNo initial amount defined for species ''%s''\n',errorMsg,speciesName);
            end           
            % Convert initialAmount to concentration and include in SBstructure
            SBstructure.parameters(numberParameters).value = initialAmount/compartmentSize;
            % Set note that type of the parameter is 'boundaey species'
            SBstructure.parameters(numberParameters).notes = '';
            % initialize type, compartment, and unittype fields
            SBstructure.parameters(numberParameters).type = 'isParameter';
            SBstructure.parameters(numberParameters).compartment = '';
            SBstructure.parameters(numberParameters).unittype = '';
            % Increment
            numberParameters = numberParameters+1;
        elseif numberRulesSpecies == 1,
            % One rule has been detected for this species
            % Boundary species may appear as products and reactants in
            % reactions but they dont have to
            if strcmp(lastRuleType,'rate'),
                % Boundary species added as a state, since rule of type 'rate'
                SBstructure.states(numberStates).name = speciesName;
                % check if an initial amount is set for a species
                if length(initialAmount)==0,
                    errorMsg = sprintf('%s\nNo initial amount defined for species ''%s''\n',errorMsg,speciesName);
                end           
                % Convert initialAmount to concentration and include in SBstructure
                SBstructure.states(numberStates).initialCondition = initialAmount/compartmentSize;
                % Include the RHS formula as the ODE for this species state
                SBstructure.states(numberStates).ODE = lastRuleFormula;
                % Set note that type of the state is 'boundary species'
                SBstructure.states(numberStates).notes = '';
                % initialize type, compartment, and unittype fields
                SBstructure.states(numberStates).type = 'isParameter';
                SBstructure.states(numberStates).compartment = '';
                SBstructure.states(numberStates).unittype = '';
                % Increment
                numberStates = numberStates+1;
            else
                % Boundary species added as a variable, since rule of type 'scalar'
                % Determine the index where to include the variable (based
                % on the ordering of the scalar rules)
                indexVariable = find(orderScalarRules==lastRuleIndex);
                SBstructure.variables(indexVariable).name = speciesName;
                SBstructure.variables(indexVariable).formula = lastRuleFormula;
                SBstructure.variables(indexVariable).notes = '';
                % initialize type, compartment, and unittype fields
                SBstructure.variables(indexVariable).type = 'isParameter';
                SBstructure.variables(indexVariable).compartment = '';
                SBstructure.variables(indexVariable).unittype = '';
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS (INCLUDE IN STRUCTURE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First process the global parameters and check if rules exist for them
%   - rate rule: as state
%   - scalar rule: as variable
%   - no rule: as parameter
%
% Then process the local parameters. They can only be parameters, no 
% rules can exist for them. In case of collisions:
%
% parameter-parameter: colliding parameters are prefixed by their reaction 
%                      names (global parameters are not prefixed)
%
% parameter-species: are prefixed by reaction name
%
% parameter-compartment: are prefixed by reaction name
%
% Start with global parameters
for k = 1:length(SBMLmodel.parameter),
    parameterName = SBMLmodel.parameter(k).name;
    parameterValue = SBMLmodel.parameter(k).value;
    % Check the number of rules the parameter is in and get the index
    % of the last rule the parameter is used in together with its type
    % ('rate' or 'scalar') and the right hand side expression 'formula'
    [numberRulesParameter,lastRuleIndex,lastRuleType,lastRuleFormula,errorMsg] = getRules(parameterName,SBMLmodel,errorMsg);
    if numberRulesParameter > 1,
        % Parameter defined by more than one rule => error
        errorMsg = sprintf('%s\nParameter ''%s'' defined by more than one rule\n',errorMsg,parameterName);
    end
    if numberRulesParameter == 0,
        % Include parameter as parameter
        SBstructure.parameters(numberParameters).name = parameterName;
        % Check if a value is set for the parameter
        if length(parameterValue)==0,
            errorMsg = sprintf('%s\nNo value defined defined for parameter ''%s''\n',errorMsg,parameterName);
        end
        % Include value in SBstructure
        SBstructure.parameters(numberParameters).value = parameterValue;
        % Set note that type of the parameter is 'global'
        SBstructure.parameters(numberParameters).notes = '';
        % initialize type, compartment, and unittype fields
        SBstructure.parameters(numberParameters).type = 'isParameter';
        SBstructure.parameters(numberParameters).compartment = '';
        SBstructure.parameters(numberParameters).unittype = '';                
        % Increment
        numberParameters = numberParameters+1;
    elseif numberRulesParameter == 1,
        % One rule has been detected for this parameter
        if strcmp(lastRuleType,'rate'),
            % Parameter added as a state, since rule of type 'rate'
            SBstructure.states(numberStates).name = parameterName;
            % check if an initial amount is set for the parameter
            if length(parameterValue)==0,
                errorMsg = sprintf('%s\nNo value defined for parameter ''%s''\n',errorMsg,parameterName);
            end
            % Include value as initial value in SBstructure
            SBstructure.states(numberStates).initialCondition = parameterValue;
            % Include the RHS formula as the ODE for this parameter state
            SBstructure.states(numberStates).ODE = lastRuleFormula;
            % Set note that type of the state is 'parameter rule'
            SBstructure.states(numberStates).notes = '';
            % initialize type, compartment, and unittype fields
            SBstructure.states(numberStates).type = 'isParameter';
            SBstructure.states(numberStates).compartment = '';
            SBstructure.states(numberStates).unittype = '';                
            % Increment
            numberStates = numberStates+1;
        else
            % Parameter added as a variable, since rule of type 'scalar'
            % Determine the index where to include the variable (based
            % on the ordering of the scalar rules)
            indexVariable = find(orderScalarRules==lastRuleIndex);
            SBstructure.variables(indexVariable).name = parameterName;
            SBstructure.variables(indexVariable).formula = lastRuleFormula;
            SBstructure.variables(indexVariable).notes = '';
            % initialize type, compartment, and unittype fields
            SBstructure.variables(indexVariable).type = 'isParameter';
            SBstructure.variables(indexVariable).compartment = '';
            SBstructure.variables(indexVariable).unittype = '';               
        end
    end
end
% Then include the local parameters within the different reactions
% If a parameter name collision has been detected, prefix the parameter
% with the reaction name it is used in 
for k1 = 1:length(SBMLmodel.reaction),
    reactionName = SBMLmodel.reaction(k1).name;
    if ~isempty(SBMLmodel.reaction(k1).kineticLaw),
        for k2 = 1:length(SBMLmodel.reaction(k1).kineticLaw.parameter),
            parameterName = SBMLmodel.reaction(k1).kineticLaw.parameter(k2).name;
            parameterValue = SBMLmodel.reaction(k1).kineticLaw.parameter(k2).value;
            % check if a value is set for the parameter
            if length(parameterValue)==0,
                errorMsg = sprintf('%s\nNo value defined for parameter ''%s'' in reaction ''%s''\n',errorMsg,parameterName,reactionName);
            end
            % add the reaction name to the parameter and update the kinetic rate law
            parameterName = char([double(reactionName) double('_') double(parameterName)]);  % fastest strcat ;)
            % Include parameter in the structure
            SBstructure.parameters(numberParameters).name = parameterName;
            SBstructure.parameters(numberParameters).value = parameterValue;
            SBstructure.parameters(numberParameters).notes = '';
            % initialize type, compartment, and unittype fields
            SBstructure.parameters(numberParameters).type = 'isParameter';
            SBstructure.parameters(numberParameters).compartment = '';
            SBstructure.parameters(numberParameters).unittype = '';
            numberParameters = numberParameters+1;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COMPARTMENTS (INCLUDE IN STRUCTURE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Include compartments
% No rule => as parameters
%
% Scalar rule => as variables
%       
% Rate rule => as states
%
for k = 1:length(SBMLmodel.compartment),
    compartmentName = SBMLmodel.compartment(k).name;
    compartmentValue = SBMLmodel.compartment(k).volume;
    % Check the number of rules the compartment is in and get the index
    % of the last rule the compartment is used in together with its type
    % ('rate' or 'scalar') and the right hand side expression 'formula'
    [numberRulesCompartment,lastRuleIndex,lastRuleType,lastRuleFormula,errorMsg] = getRules(compartmentName,SBMLmodel,errorMsg);
    if numberRulesCompartment > 1,
        % Compartment defined by more than one rule => error
        errorMsg = sprintf('%s\nCompartment ''%s'' defined by more than one rule\n',errorMsg,compartmentName);
    end
    if numberRulesCompartment == 0,
        % Include compartment as parameter
        SBstructure.parameters(numberParameters).name = compartmentName;
        % Check if a value is set for the parameter
        if length(compartmentValue)==0,
            errorMsg = sprintf('%s\nNo value defined defined for compartment ''%s''\n',errorMsg,compartmentName);
        end
        % Include value in SBstructure
        SBstructure.parameters(numberParameters).value = compartmentValue;
        % Set note that type of the parameter is 'compartment size'
        SBstructure.parameters(numberParameters).notes = '';
        % initialize type, compartment, and unittype fields
        SBstructure.parameters(numberParameters).type = 'isCompartment';
        SBstructure.parameters(numberParameters).compartment = SBMLmodel.compartment(k).outside;
        SBstructure.parameters(numberParameters).unittype = '';
        % Increment
        numberParameters = numberParameters+1;
    elseif numberRulesCompartment == 1,
        % One rule has been detected for this compartment
        if strcmp(lastRuleType,'rate'),
            % Compartment added as a state, since rule of type 'rate'
            SBstructure.states(numberStates).name = compartmentName;
            % check if an initial amount is set for the compartment
            if length(compartmentValue)==0,
                errorMsg = sprintf('%s\nNo value defined for compartment ''%s''\n',errorMsg,compartmentName);
            end
            % Include value as initial value in SBstructure
            SBstructure.states(numberStates).initialCondition = compartmentValue;
            % Include the RHS formula as the ODE for this compartment state
            SBstructure.states(numberStates).ODE = lastRuleFormula;
            % Set note that type of the state is 'compartment size'
            SBstructure.states(numberStates).notes = '';
            % initialize type, compartment, and unittype fields
            SBstructure.states(numberStates).type = 'isCompartment';
            SBstructure.states(numberStates).compartment = SBMLmodel.compartment(k).outside;
            SBstructure.states(numberStates).unittype = '';
            % Increment
            numberStates = numberStates+1;
            % compartment_rate variable right hand side
            compartmentDotRHS = lastRuleFormula;
        else
            % Compartment added as a variable, since rule of type 'scalar'
            % Determine the index where to include the variable (based
            % on the ordering of the scalar rules)
            indexVariable = find(orderScalarRules==lastRuleIndex);
            SBstructure.variables(indexVariable).name = compartmentName;
            SBstructure.variables(indexVariable).formula = lastRuleFormula;
            SBstructure.variables(indexVariable).notes = '';
            % initialize type, compartment, and unittype fields
            SBstructure.variables(indexVariable).type = 'isCompartment';
            SBstructure.variables(indexVariable).compartment = SBMLmodel.compartment(k).outside;
            SBstructure.variables(indexVariable).unittype = '';
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REACTIONS (INCLUDE IN STRUCTURE)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read in reaction information into the model structure. Replace colliding 
% parameter names with their prefixed versions (reaction name). Reactions 
% are included into the model structure in the 'reactions' field. This leads 
% to shorter ODE definitions for the specie states, which then contain the 
% reaction names with the correct stoichiometries. The 'reactions' field
% contains three fields: 'name', 'formula', 'notes'
%
% Cycle through all reactions
for k1 = 1:length(SBMLmodel.reaction),
    reactionName = SBMLmodel.reaction(k1).name;
    kineticLaw = SBMLmodel.reaction(k1).kineticLaw;
    % check if a kineticLaw is given
    if length(kineticLaw)==0,
        errorMsg = sprintf('%s\nNo kinetic law is given for reaction ''%s''\n',errorMsg,reactionName);
    end
    % Process the kineticLaw formula by replacing the local parameter names
    % by reactionname_parameternames and replace the mathml expressions
    if ~isempty(kineticLaw),
        formula = kineticLaw.formula;
        if length(kineticLaw.parameter) > 0,
            for k2 = 1:length(kineticLaw.parameter),
                parameterName = kineticLaw.parameter(k2).name;
                newParameterName = char([double(reactionName) double('_') double(parameterName)]); % fastest strcat
                formula = exchangeStringInString(formula,parameterName,newParameterName);
            end
        end
    else
        formula = '';
    end
    % replace MathML expressions in formula against MATLAB expressions
    [formula] = replaceMathMLexpressions(formula);
    % Include reaction information in the structure
    SBstructure.reactions(k1).name = reactionName;
    SBstructure.reactions(k1).formula = formula;
    SBstructure.reactions(k1).notes = '';
    SBstructure.reactions(k1).reversible = SBMLmodel.reaction(k1).reversible;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS REACTION ODES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add the reaction variable names with correct stoichiometries to the 'ODE' field
% for the correct species. The 'correct' species are only non-boundary and 
% non-constant species for which no rules exist. These are listed in
% the 'speciesODElist'
%
for k1 = 1:length(speciesODElist),
    % get name of species and the index of the corresponding state
    speciesName = speciesODElist(k1).name;
    stateIndex = speciesODElist(k1).stateIndex;
    % cycle through all the reactions and look in which reactions the
    % species is changed
    for k2 = 1:length(SBstructure.reactions),
        % get reaction name for inclusion in ODEstring
        reactionName = SBstructure.reactions(k2).name;
        % cycle through the reactants of the current reaction
        % (ordering of reactions in SBstructure and SBMLmodel are equal)
        for k3 = 1:length(SBMLmodel.reaction(k2).reactant),
            reactantName = SBMLmodel.reaction(k2).reactant(k3).species;
            reactantStoichiometry = abs(double(SBMLmodel.reaction(k2).reactant(k3).stoichiometry) / double(SBMLmodel.reaction(k2).reactant(k3).denominator));
            % If species name and reactant name are equal then add reaction
            % to species ODE
            if strcmp(reactantName,speciesName),
                % construct the string to add to the ODE of the current reactant
                if reactantStoichiometry > 1,
                    ODEstringAdd = sprintf('-%s*%s',num2str(reactantStoichiometry),reactionName);
                else
                    ODEstringAdd = sprintf('-%s',reactionName);
                end
                % add the ODEstringAdd to the 'ODE' field of this species
                SBstructure.states(stateIndex).ODE = strcat(SBstructure.states(stateIndex).ODE,ODEstringAdd);
                break;
            end
        end
        % cycle through the products of the current reaction
        % (ordering of reactions in SBstructure and SBMLmodel are equal)
        for k3 = 1:length(SBMLmodel.reaction(k2).product),
            productName = SBMLmodel.reaction(k2).product(k3).species;
            productStoichiometry = abs(double(SBMLmodel.reaction(k2).product(k3).stoichiometry) / double(SBMLmodel.reaction(k2).product(k3).denominator));
            % If species name and product name are equal then add reaction
            % to species ODE
            if strcmp(productName,speciesName),
                % construct the string to add to the ODE of the current product
                if productStoichiometry > 1,
                    ODEstringAdd = sprintf('+%s*%s',num2str(productStoichiometry),reactionName);
                else
                    ODEstringAdd = sprintf('+%s',reactionName);
                end
                % add the ODEstringAdd to the 'ODE' field of this species
                SBstructure.states(stateIndex).ODE = strcat(SBstructure.states(stateIndex).ODE,ODEstringAdd);
                break;
            end
        end
    end
end
% Now all ODEs should be defined - check if everything is correct!
% Finally cycle through all the states and check if all ODEs have been defined
for k = 1:length(SBstructure.states),
    stateName = SBstructure.states(k).name;
    stateODE = SBstructure.states(k).ODE;
    if length(stateODE) == 0,
        % Issue a warning only and set ODE to 0
        disp(sprintf('No ODE defined for state ''%s''',stateName));
        SBstructure.states(k).ODE = '0';
    end
end
% After initial construction of the ODEs their units need to be changed.
% Since species are in concentrations and the rate laws from the reactions
% return amount/time it is necessary to convert the odes from amount/time to
% concentration/time. This is done here.
%
% We differentiate the case where:
%   - single compartment model with constant compartment size = 1
%     (nothing is done additionally)
%   - multi compartment model or single compartment model with compartment 
%     size different from one. In this case each complete ode is divided by the
%     compartment size the species is in
%
for k1 = 1:length(speciesODElist),
    % get name of species and the index of the corresponding state
    speciesName = speciesODElist(k1).name;
    stateIndex = speciesODElist(k1).stateIndex;
    % check if the model has only one compartment with volume 1 and 
    % constant volume => if yes then leave ODE as it is.
    convertODE = 1;
    if length(SBMLmodel.compartment) == 1,
        % only one compartment
        if SBMLmodel.compartment(1).volume == 1,
            % size of compartment is 1
            for k = 1:length(SBstructure.parameters),
                % check is compartment size is constant (then it is defined as
                % a parameter)
                if strcmp(SBMLmodel.compartment(1).name,SBstructure.parameters(k).name),
                    % compartment defined as a parameter and thus no species
                    % names need to be exchanged
                    convertODE = 0;
                    % return directly to the main function
                end
            end
        end
    end
    if convertODE,
        % convert the ODE expression
        % get compartment name for species
        compartmentName = getCompartmentNameSpecies(speciesName,SBMLmodel);
        % construct new ODE string
        ODEstring = SBstructure.states(stateIndex).ODE;
        newODEstring = sprintf('(%s)/%s',ODEstring,compartmentName);
        % include the new ODE string in structure
        SBstructure.states(stateIndex).ODE = newODEstring;
    end
end
% return from main function
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIND RULES FOR A CERTAIN VARIABLE NAME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cycle through all the rules and check if the given variable name
% (species or parameter or compartment) appears as left hand side (LHS) in a rule
%
% numberRulesVariable: number of rules the name appears in as LHS
% lastRuleIndex: the index of the rule the name was last detected
% lastRuleType: 'scalar' or 'rate'
% lastRuleFormula: the RHS expression for the variable
%
% We assume that the name of the variable can appear in the
% 'variable' field, in the 'species' field, in the 'compartment' field,
% or in the 'name' field of the 'rule' field. 
% Some SBML example models had species in the 'variable' field!
%
function [numberRulesVariable,lastRuleIndex,lastRuleType,lastRuleFormula,errorMsg] = getRules(variableName,SBMLmodel,errorMsg)
% Initialize return values
numberRulesVariable = 0;
lastRuleIndex = [];
lastRuleType = '';
lastRuleFormula = '';
for k = 1:length(SBMLmodel.rule),
    if strcmp(variableName,SBMLmodel.rule(k).variable) || strcmp(variableName,SBMLmodel.rule(k).species) || strcmp(variableName,SBMLmodel.rule(k).name) || strcmp(variableName,SBMLmodel.rule(k).compartment),
        % The name was found as LHS in a rule
        numberRulesVariable = numberRulesVariable+1;
        lastRuleIndex = k;
        if ~isempty(strfind(SBMLmodel.rule(k).typecode,'RATE')) || ~isempty(strfind(SBMLmodel.rule(k).type,'rate')),
            % Rule is a rate rule
            lastRuleType = 'rate';
        else
            % Rule is a scalar rule
            lastRuleType = 'scalar';
        end
        % get the RHS expression of the rule
        lastRuleFormula = SBMLmodel.rule(k).formula;
    end
end
% replace MathML expressions against MATLAB expressions
[lastRuleFormula] = replaceMathMLexpressions(lastRuleFormula);
%return
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK IF A CERTAIN SPECIES APPEARS AS PRODUCT AND/OR REACTANT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% In the case that a non-boundary species is defined by a rule
% we need to check if it is also present as product or reactant in any
% reaction. In this case the SBML model is not correctly defined and
% reactionSpeciesPresent = 1
function [reactionSpeciesPresent] = checkReactionsForSpecies(speciesName,SBMLmodel)
% Initialize return value
reactionSpeciesPresent = 0;
reactionComponents = {};
for k = 1:length(SBMLmodel.reaction),
    reactionComponents = {reactionComponents{:},SBMLmodel.reaction(k).reactant.species,SBMLmodel.reaction(k).product.species};
end
if ~isempty(strmatch(speciesName,reactionComponents,'exact')),
    reactionSpeciesPresent = 1;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DELETE WHITESPACES IN STRINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Useful for taking away whitespaces in kineticLaw formulas, as
% seen in some example models
function [outputString] = removeWhiteSpace(inputString)
outputString = strrep(inputString,' ','');
% return
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET COMPARTMENT NAME FOR SPECIES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [compartmentName] = getCompartmentNameSpecies(species,SBMLmodel)
compartmentName = SBMLmodel.species(strmatch(species,{SBMLmodel.species.name},'exact')).compartment;
% return
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXCHANGE STRING IN STRING - WITH CHECK OF ALLOWED CHARACTERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find occurrences of stringOld in fullString and exchange it with
% stringNew. 
function [processedString] = exchangeStringInString(fullString,stringOld,stringNew);
% really simple!!!
processedString = regexprep(fullString,strcat('\<',stringOld,'\>'),stringNew);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXCHANGE MATHML FUNCTION NAMES IN ALL FORMULAS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MathML functions can have different names than MATLAB function names
% we need to exchange them in the formula strings. The term 'MathML
% expression' relates to what is returned by TranslateSBML!
function [newFormula] = replaceMathMLexpressions(formula)
% MathML expressions that need simple exchange with corresponding MATLAB
% expressions:
MathMLexpressions = {'\<piecewise\>','\<and\>','\<or\>','\<arccos\>','\<arcsin\>','\<arctan\>','\<ceiling\>','\<ln\>',...
    '\<pow\>','\<arccos\>','\<arccosh\>','\<arccot\>','\<arccoth\>','\<arccsc\>','\<arccsch\>',...
    '\<arcsec\>','\<arcsech\>','\<arcsin\>','\<arcsinh\>','\<arctan\>','\<arctanh\>',...
    '\<exponentiale'};
MATLABexpressions = {'piecewiseSB','andSB','orSB','acos','asin','atan','ceil','log', ...
    'power','acos','acosh','acot','acoth','acsc','acsch',...
    'asec','asech','asin','asinh','atan','atanh'...
    'exp(1)'};
% find indices
newFormula = regexprep(formula,MathMLexpressions,MATLABexpressions);
return








