function [namesOK] = checkNamesInSBmodel(model, varargin)
% checkNamesInSBmodel
% checks the used names in the given SBmodel (empty name fields and names
% used more than once)
%
%
% USAGE:
% ======
% [namesOK] = checkSBmodelNamesForDoublets(model)
%
% model: SBmodel where names are to test
% 
% namesOK: true if no errors in name definitions where found
%          otherwise false
%

% Information:
% ============
% Systems Biology Toolbox for MATLAB
% Copyright (C) 2005-2007 Fraunhofer-Chalmers Research Centre
% Author of the toolbox: Henning Schmidt, henning@hschmidt.de
% Programmed by Gunnar Drews, Student at University of Rostock
% Department of Computerscience, gunnar.drews@uni-rostock.de
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

silentFlag = 0;
if (nargin == 2)
    silentFlag = varargin{1};
end


namesOK = false;
emptyNames = {};
emptyNamesIndex = 1;
nameError = false;
errorGroup = {};
errorName = '';
nameDoubletError = false;
modelStruct = SBstruct(model);
variableCount = length(modelStruct.variables);

% check names of SBmodel states 
stateNames = {};
index = 1;
for n = 1 : length(modelStruct.states),
    name = modelStruct.states(n).name;
    if strcmp(name, ''),
        emptyNames{emptyNamesIndex} = char([double('State No. '), double(numstr(n))]);
        emptyNamesIndex = emptyNamesIndex + 1;
        nameError = true;
    else
        if matchStringOnArray(name, stateNames),
            errorGroup{1} = 'States';
            errorGroup{2} = 'States';
            errorName = name;
            nameDoubletError = true;
        else
            stateNames{index} = name;
            index = index + 1;
        end
    end
end
% check names of SBmodel parameters
parameterNames = {};
index = 1;
if ~nameDoubletError,
    for n = 1 : length(modelStruct.parameters),
        name = modelStruct.parameters(n).name;
        if strcmp(name, ''),
            emptyNames{emptyNamesIndex} = char([double('Parameter No. '), double(numstr(n))]);
            emptyNamesIndex = emptyNamesIndex + 1;
            nameError = true;
        else
            if matchStringOnArray(name, stateNames),
                errorGroup{1} = 'States';
                errorGroup{2} = 'Parameters';
                errorName = name;
                nameDoubletError = true;
            elseif matchStringOnArray(name, parameterNames),
                errorGroup{1} = 'Parameters';
                errorGroup{2} = 'Parameters';
                errorName = name;
                nameDoubletError = true;
            else
                parameterNames{index} = name;
                index = index + 1;
            end
        end
    end
end
% check names of SBmodel variables
variableNames = {};
index = 1;
if ~nameDoubletError,
    for n = 1 : length(modelStruct.variables),
        name = modelStruct.variables(n).name;
        if strcmp(name, ''),
            emptyNames{emptyNamesIndex} = char([double('Variable No. '), double(numstr(n))]);
            emptyNamesIndex = emptyNamesIndex + 1;
            nameError = true;
        else
            if matchStringOnArray(name, stateNames),
                errorGroup{1} = 'States';
                errorGroup{2} = 'Variables';
                errorName = name;
                nameDoubletError = true;
            elseif matchStringOnArray(name, parameterNames),
                errorGroup{1} = 'Parameters';
                errorGroup{2} = 'Variables';
                errorName = name;
                nameDoubletError = true;
            elseif matchStringOnArray(name, variableNames),
                errorGroup{1} = 'Variables';
                errorGroup{2} = 'Variables';
                errorName = name;
                nameDoubletError = true;
            else
                variableNames{index} = name;
                index = index + 1;
            end
        end
    end
end
% check names of SBmodel reactions
reactionNames = {};
index = 1;
if ~nameDoubletError,
    for n = 1 : length(modelStruct.reactions),
        name = modelStruct.reactions(n).name;
        if strcmp(name, ''),
            emptyNames{emptyNamesIndex} = char([double('Reaction No. '), double(numstr(n))]);
            emptyNamesIndex = emptyNamesIndex + 1;
            nameError = true;
        else
            if matchStringOnArray(name, stateNames),
                errorGroup{1} = 'States';
                errorGroup{2} = 'Reactions';
                errorName = name;
                nameDoubletError = true;
            elseif matchStringOnArray(name, parameterNames),
                errorGroup{1} = 'Parameters';
                errorGroup{2} = 'Reactions';
                errorName = name;
                nameDoubletError = true;
            elseif matchStringOnArray(name, variableNames),
                errorGroup{1} = 'Variables';
                errorGroup{2} = 'Reactions';
                errorName = name;
                nameDoubletError = true;
            elseif matchStringOnArray(name, reactionNames),
                errorGroup{1} = 'Reactions';
                errorGroup{2} = 'Reactions';
                errorName = name;
                nameDoubletError = true;
            else
                reactionNames{index} = name;
                index = index + 1;
            end
        end
    end
end
% check names of SBmodel functions
functionNames = {};
index = 1;
if ~nameDoubletError,
    for n = 1 : length(modelStruct.functions),
        name = modelStruct.functions(n).name;
        if strcmp(name, ''),
            emptyNames{emptyNamesIndex} = char([double('Function No. '), double(numstr(n))]);
            emptyNamesIndex = emptyNamesIndex + 1;
            nameError = true;
        else
            if matchStringOnArray(name, stateNames),
                errorGroup{1} = 'States';
                errorGroup{2} = 'Functions';
                errorName = name;
                nameDoubletError = true;
            elseif matchStringOnArray(name, parameterNames),
                errorGroup{1} = 'Parameters';
                errorGroup{2} = 'Functions';
                errorName = name;
                nameDoubletError = true;
            elseif matchStringOnArray(name, variableNames),
                errorGroup{1} = 'Variables';
                errorGroup{2} = 'Functions';
                errorName = name;
                nameDoubletError = true;
            elseif matchStringOnArray(name, reactionNames),
                errorGroup{1} = 'Reactions';
                errorGroup{2} = 'Functions';
                errorName = name;
                nameDoubletError = true;
            elseif matchStringOnArray(name, functionNames),
                errorGroup{1} = 'Reactions';
                errorGroup{2} = 'Functions';
                errorName = name;
                nameDoubletError = true;
            else
                functionNames{index} = name;
                index = index + 1;
            end
        end
    end
end
% check names of SBmodel events
eventNames = {};
index = 1;
if ~nameDoubletError,
    for n = 1 : length(modelStruct.events),
        name = modelStruct.events(n).name;
        if strcmp(name, ''),
            emptyNames{emptyNamesIndex} = char([double('Event No. '), double(numstr(n))]);
            emptyNamesIndex = emptyNamesIndex + 1;
            nameError = true;
        else
            if matchStringOnArray(name, stateNames),
                errorGroup{1} = 'States';
                errorGroup{2} = 'Events';
                errorName = name;
                nameDoubletError = true;
            elseif matchStringOnArray(name, parameterNames),
                errorGroup{1} = 'Parameters';
                errorGroup{2} = 'Events';
                errorName = name;
                nameDoubletError = true;
            elseif matchStringOnArray(name, variableNames),
                errorGroup{1} = 'Variables';
                errorGroup{2} = 'Events';
                errorName = name;
                nameDoubletError = true;
            elseif matchStringOnArray(name, reactionNames),
                errorGroup{1} = 'Reactions';
                errorGroup{2} = 'Events';
                errorName = name;
                nameDoubletError = true;
            elseif matchStringOnArray(name, functionNames),
                errorGroup{1} = 'Reactions';
                errorGroup{2} = 'Events';
                errorName = name;
                nameDoubletError = true;
            elseif matchStringOnArray(name, eventNames),
                errorGroup{1} = 'Events';
                errorGroup{2} = 'Events';
                errorName = name;
                nameDoubletError = true;
            else
                eventNames{index} = name;
                index = index + 1;
            end
        end
    end
end

errorMessage = {};
errorIndex = 1;
if nameError,
    errorMessage{errorIndex} = 'The following SBmodel components have no defined name(s):';
    errorIndex = errorIndex + 1;
    for n = 1 : length(emptyNames),
        errorMessage{errorIndex} = emptyNames(n);
        errorIndex = errorIndex + 1;
    end
    errorMessage{errorIndex} = '';
    errorIndex = errorIndex + 1;
end
if nameDoubletError,
    errorMessage{errorIndex} = char([double('The SBmodel component "'), double(errorName), double('" has been found ')]);
    errorIndex = errorIndex + 1;
    if strcmp(errorGroup{1}, errorGroup{2}),
        errorMessage{errorIndex} = char([double(' twice in SBmodel '), double(errorGroup{1}), double('.')]);
        errorIndex = errorIndex + 1;
    else
        errorMessage{errorIndex} = char([double(' in SBmodel '), double(errorGroup{1}), double('and SBmodel '), double(errorGroup{2}), double('.')]);
        errorIndex = errorIndex + 1;
    end
    errorMessage{errorIndex} = '';
    errorIndex = errorIndex + 1;
end

if (nameError || nameDoubletError),
    if ~silentFlag,
        messageOutput(errorMessage, 1);
    end
else    
    namesOK = true;
end

return