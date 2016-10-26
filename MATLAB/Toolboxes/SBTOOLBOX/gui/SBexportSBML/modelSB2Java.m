function [SBmodelJava] = modelSB2Java(model)
% modelSB2Java 
%
% converts a SBmodel to a SBmodelJava object which is used for the
% SBexportSBML gui.
%
% USAGE:
% ======
% [SBmodelJava] = modelSB2Java(model)
%
% model: SBmodel to convert to a Java object

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

% get the MATLAB structure of the provided model
sbm = struct(model);
% create an empty SBmodelJava object
sbmj = auxiliary.javamodel.SBmodelJava;

% copy SBmodel name and notes to the Java object
%if exist(sbm.name),
    sbmj.setName(sbm.name);
%end
%if exist(sbm.notes),
    sbmj.setNotes(sbm.notes);
%end

% copy all functions from the SBmodel to the Java object
%if exist('sbm.functions', 'var'),
    functionsCount=length(sbm.functions);
    if functionsCount>0,
        % we need to set the number of functions before we can access
        % the specific fields of the object
        sbmj.setNumberOfFunctions(functionsCount);
        for k = 1 : functionsCount,
            % count variable "k" has to be decreased by one because Java
            % arrays start with zero
            % struct2cell array conversion seems to be the easiest way for information
            % passing because Java can handle MATLAB cell arrays or autocasting
            % for cell arrays within MATLAB works good enough
            %if (exist('sbm.functions(k).name', 'var') && exist('sbm.functions(k).arguments', 'var') && exist('sbm.functions(k).formula', 'var') && exist('sbm.functions(k).notes', 'var')),
                sbmj.setFunction(k-1, struct2cell(sbm.functions(k)));
            %end
        end
    end
%end

% copy all states from the SBmodel to the Java object
%if exist('sbm.states', 'var'),
    statesCount=length(sbm.states);
    if statesCount>0,
        % we need to set the number of states before we can access
        % the specific fields of the object
        sbmj.setNumberOfStates(statesCount);
        for k = 1 : statesCount,
            % count variable "k" has to be decreased by one because Java
            % arrays start with zero
            % struct2cell array conversion seems to be the easiest way for information
            % passing because Java can handle MATLAB cell arrays or autocasting
            % for cell arrays within MATLAB works good enough
            %if exist('sbm.states(k).name', 'var') && exist('sbm.states(k).initialCondition', 'var') && exist('sbm.states(k).ODE', 'var') && exist('sbm.states(k).type', 'var') && exist('sbm.states(k).compartment', 'var') && exist('sbm.states(k).unittype', 'var') && exist('sbm.states(k).notes', 'var'),
                sbmj.setState(k-1, struct2cell(sbm.states(k)));
            %end
        end
    end
%end

% copy all parameters from the SBmodel to the Java object
%if exist('sbm.parameters', 'var'),
    parametersCount=length(sbm.parameters);
    if parametersCount>0,
        % we need to set the number of parameters before we can access
        % the specific fields
        sbmj.setNumberOfParameters(parametersCount);
        for k = 1 : parametersCount,
            % count variable "k" has to be decreased by one because Java
            % arrays start with zero
            % struct2cell array conversion seems to be the easiest way for information
            % passing because Java can handle MATLAB cell arrays or autocasting
            % for cell arrays within MATLAB works good enough
           %if exist('sbm.parameters(k).name', 'var') && exist('sbm.parameters(k).value', 'var') && exist('sbm.parameters(k).type', 'var') && exist('sbm.parameters(k).compartment', 'var') && exist('sbm.parameters(k).unittype', 'var') && exist('sbm.parameters(k).notes', 'var'),
                sbmj.setParameter(k-1, struct2cell(sbm.parameters(k)));
            %end
        end
    end
%end

% copy all variables from the SBmodel to the Java object
%if exist('sbm.variables', 'var'),
    variablesCount=length(sbm.variables);
    if variablesCount>0,
        % we need to set the number of variables before we can access
        % the specific fields of the object
        sbmj.setNumberOfVariables(variablesCount);
        for k = 1 : variablesCount,
            % count variable "k" has to be decreased by one because Java
            % arrays start with zero
            % struct2cell array conversion seems to be the easiest way for information
            % passing because Java can handle MATLAB cell arrays or autocasting
            % for cell arrays within MATLAB works good enough
            %if exist('sbm.variables(k).name', 'var') && exist('sbm.variables(k).formula', 'var') && exist('sbm.variables(k).type', 'var') && exist('sbm.variables(k).compartment', 'var') && exist('sbm.variables(k).unittype', 'var') && exist('sbm.variables(k).notes', 'var'),
                sbmj.setVariable(k-1, struct2cell(sbm.variables(k)));
            %end
        end
    end
%end

% copy all reactions from the SBmodel to the Java object
%if exist('sbm.reactions', 'var'),
    reactionsCount=length(sbm.reactions);
    if reactionsCount>0,
        % we need to set the number of reactions before we can access
        % the specific fields of the object
        sbmj.setNumberOfReactions(reactionsCount);
        for k = 1 : reactionsCount,
            % count variable "k" has to be decreased by one because Java
            % arrays start with zero
            % struct2cell array conversion seems to be the easiest way for information
            % passing because Java can handle MATLAB cell arrays or autocasting
            % for cell arrays within MATLAB works good enough
            %if exist('sbm.reactions(k).name', 'var') && exist('sbm.reactions(k).formula', 'var') && exist('sbm.reactions(k).notes', 'var') && exist('sbm.reactions(k).reversible', 'var'),
                sbmj.setReaction(k-1, struct2cell(sbm.reactions(k)));
            %end
        end
    end
%end

% copy all events from the SBmodel to the Java object
%if exist('sbm.events', 'var'),
    eventsCount=length(sbm.events);
    if eventsCount>0,
        % we need to set the number of events before we can access
        % the specific fields of the object
        sbmj.setNumberOfEvents(eventsCount);
        for k = 1 : eventsCount,
            % count variable "k" has to be decreased by one because Java
            % arrays start with zero
            % struct2cell array conversion seems to be the easiest way for information
            % passing because Java can handle MATLAB cell arrays or autocasting
            % for cell arrays within MATLAB works good enough
            %if exist('sbm.events(k).name', 'var'),
                sbmj.getEvent(k-1).setName(sbm.events(k).name);
            %end
            %if exist('sbm.events(k).trigger', 'var'),
                sbmj.getEvent(k-1).setTrigger(sbm.events(k).trigger);
            %end
            %if exist('sbm.events(k).assignment', 'var'),
                assignmentCount=length(sbm.events(k).assignment);
                if assignmentCount>0,
                    % events contain another array (assignment)
                    % we need to set the number of events before we can access
                    % the specific fields of the object
                    sbmj.getEvent(k-1).setNumberOfAssignments(assignmentCount);
                    for k2 = 1 : assignmentCount,
                        %if exist('sbm.events(k).assignment(k2).variable', 'var') && exist('sbm.events(k).assignment(k2).formula', 'var'),
                            sbmj.getEvent(k-1).setAssignment(k2-1, struct2cell(sbm.events(k).assignment(k2)));
                        %end
                    end
                end
            %end
            %if exist('sbm.events(k).notes', 'var'),
                sbmj.getEvent(k-1).setNotes(sbm.events(k).notes);
            %end
        end
    end
%end

% copy FunctionsMATLAB to the Java object
%if exist('sbm.functionsMATLAB', 'var'),
    sbmj.setFunctionsMATLAB(sbm.functionsMATLAB);
%end
% take local Java object sbmj and set it as return parameter
SBmodelJava=sbmj;
return