function [sbm] = java2modelSB(model)
% java2modelSB 
%
% converts a SBmodelJava object to the usual SBmodel.
%
% USAGE:
% ======
% [SBmodel] = java2modelSB(model)
%
% model: SBmodelJava object to convert to the usual SBmodel

% Information:
% ============
% Systems Biology Toolbox for MATLAB
% Copyright (C) 2005-2007 Fraunhofer-Chalmers Research Centre
% Author of the toolbox: Henning Schmidt, henning@hschmidt.de
% Programmed by Gunnar Drews, gunnar.drews@uni-rostock.de
% Student at University of Rostock Department of Computerscience
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

% copy SBmodelJava name and notes to the SBmodel structure
SBstructure.name = char(model.getName());
SBstructure.notes = char(model.getNotes());

% copy all functions from the SBmodelJava to the SBmodel structure
functionsCount=model.getNumberOfFunctions();
if functionsCount>0,
    for k = 1 : functionsCount,
        % for the Java object the count variable "k" has to be decreased by
        % one because Java arrays start with zero
        SBstructure.functions(k).name = char(model.getFunction(k-1).getName());
        SBstructure.functions(k).arguments = char(model.getFunction(k-1).getArguments());
        SBstructure.functions(k).formula = char(model.getFunction(k-1).getFormula());
        SBstructure.functions(k).notes = char(model.getFunction(k-1).getNotes());
    end
end

% copy all states from the SBmodelJava to the SBmodel structure
statesCount=model.getNumberOfStates();
if statesCount>0,
    for k = 1 : statesCount,
        % for the Java object the count variable "k" has to be decreased by
        % one because Java arrays start with zero
        SBstructure.states(k).name = char(model.getState(k-1).getName());
        SBstructure.states(k).initialCondition = double(model.getState(k-1).getInitialCondition());
        SBstructure.states(k).ODE = char(model.getState(k-1).getODE());
        SBstructure.states(k).type = char(model.getState(k-1).getType().asSBToolboxExpression());
        SBstructure.states(k).compartment = char(model.getState(k-1).getCompartment());
        SBstructure.states(k).unittype = char(model.getState(k-1).getUnitType().asSBToolboxExpression());
        SBstructure.states(k).notes = char(model.getState(k-1).getNotes());
    end
end

% copy all parameters from the SBmodel to the Java object
parametersCount=model.getNumberOfParameters();
if parametersCount>0,
    for k = 1 : parametersCount,
        % for the Java object the count variable "k" has to be decreased by
        % one because Java arrays start with zero
        SBstructure.parameters(k).name = char(model.getParameter(k-1).getName());
        SBstructure.parameters(k).value = double(model.getParameter(k-1).getValue());
        SBstructure.parameters(k).type = char(model.getParameter(k-1).getType().asSBToolboxExpression());
        SBstructure.parameters(k).compartment = char(model.getParameter(k-1).getCompartment());
        SBstructure.parameters(k).unittype = char(model.getParameter(k-1).getUnitType().asSBToolboxExpression());
        SBstructure.parameters(k).notes = char(model.getParameter(k-1).getNotes());
    end
end

% copy all variables from the SBmodel to the Java object
variablesCount=model.getNumberOfVariables();
if variablesCount>0,
    for k = 1 : variablesCount,
        % for the Java object the count variable "k" has to be decreased by
        % one because Java arrays start with zero
        SBstructure.variables(k).name = char(model.getVariable(k-1).getName());
        SBstructure.variables(k).formula = char(model.getVariable(k-1).getFormula());
        SBstructure.variables(k).type = char(model.getVariable(k-1).getType().asSBToolboxExpression());
        SBstructure.variables(k).compartment = char(model.getVariable(k-1).getCompartment());
        SBstructure.variables(k).unittype = char(model.getVariable(k-1).getUnitType().asSBToolboxExpression());
        SBstructure.variables(k).notes = char(model.getVariable(k-1).getNotes());
    end
end

% copy all reactions from the SBmodel to the Java object
reactionsCount=model.getNumberOfReactions();
if reactionsCount>0,
    for k = 1 : reactionsCount,
        % for the Java object the count variable "k" has to be decreased by
        % one because Java arrays start with zero
        SBstructure.reactions(k).name = char(model.getReaction(k-1).getName());
        SBstructure.reactions(k).formula = char(model.getReaction(k-1).getFormula());
        SBstructure.reactions(k).notes = char(model.getReaction(k-1).getNotes());
        SBstructure.reactions(k).reversible = double(model.getReaction(k-1).getReversible());
    end
end

% copy all events from the SBmodel to the Java object
eventsCount=model.getNumberOfEvents();
if eventsCount>0,
    for k = 1 : eventsCount,
        % for the Java object the count variable "k" has to be decreased by
        % one because Java arrays start with zero
        SBstructure.events(k).name = char(model.getEvent(k-1).getName());
        SBstructure.events(k).trigger = char(model.getEvent(k-1).getTrigger());
        assignmentCount=model.getEvent(k-1).getNumberOfAssignments();
        if assignmentCount>0,
            for k2 = 1 : assignmentCount,
                SBstructure.events(k).assignment(k2).variable = char(model.getEvent(k-1).getAssignment(k2-1).getVariable());
                SBstructure.events(k).assignment(k2).formula = char(model.getEvent(k-1).getAssignment(k2-1).getFormula());
            end
        end
        SBstructure.events(k).notes = char(model.getEvent(k-1).getNotes());
    end
end

% copy FunctionsMATLAB to the SBmodel
SBstructure.functionsMATLAB = char(model.getFunctionsMATLAB());
% take local SBstructure and set it as return parameter
sbm = SBmodel(SBstructure);
return