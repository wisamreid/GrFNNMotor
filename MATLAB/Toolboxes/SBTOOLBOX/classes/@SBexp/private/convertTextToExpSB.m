function [SBstructure,errorMsg] = convertTextToExpSB(expText)
% convertTextToExpSB: Converts a text description of an SBexp object to 
% the internal SBexp data structure representation.

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

% initialize variables
errorMsg = '';
errorConditions = '';
errorParameterChanges = '';
errorStateEvents = '';
SBstructure = [];

% cut text into pieces
expTextStructure = getPartsFromCompleteTextExpSB(expText);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SBstructure.name = strtrim(removeCharacters(expTextStructure.name));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Notes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SBstructure.notes = strtrim(expTextStructure.notes);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Conditions (parameters and initial conditions)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
try
    [SBstructure.initialconditions, SBstructure.parametersettings, errorConditions] = getConditions(expTextStructure.conditions);
catch
    errorMsg = sprintf('%sPlease check the syntax of the ''Parameter and initial conditions'' definitions.\n',errorMsg);
end
if ~isempty(errorConditions),
    errorMsg = sprintf('%s%s\n',errorMsg,errorConditions);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% Parameter changes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%try
    [SBstructure.parameterchanges, errorParameterChanges] = getParameterChanges(expTextStructure.parameterchanges);
%catch
%    errorMsg = sprintf('%sPlease check the syntax of the ''Parameter Changes'' definitions.\n',errorMsg);
%end
if ~isempty(errorParameterChanges),
    errorMsg = sprintf('%s%s\n',errorMsg,errorParameterChanges);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% State changes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%try
    [SBstructure.stateevents, errorStateEvents] = getStateEvents(expTextStructure.stateevents);
%catch
%    errorMsg = sprintf('%sPlease check the syntax of the ''States Changes'' definitions.\n',errorMsg);
%end
if ~isempty(errorStateEvents),
    errorMsg = sprintf('%s%s\n',errorMsg,errorStateEvents);
end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [initialconditions, parametersettings, error] = getConditions(conditionsString)
error = '';
% create empty structures
parametersettings = struct('name',{},'formula',{},'notes',{});
initialconditions = struct('name',{},'formula',{},'notes',{});    
% parse the strings and fill the structures
conditionsString = char([double(conditionsString) 10]); % need to add a '\n' at the end
allelements = regexp(conditionsString,'([^\n]*)=([^\n]*)\n','tokens');
for k = 1:length(allelements),
    element = allelements{k};
    leftside = strtrim(removeCharacters(element{1}));
    rightsideandcomment = strtrim(removeCharacters(element{2}));
    % parse right side and comment (remove (0) from the right side)
    rightsideandcomment = regexp(rightsideandcomment,'([^%]*)','tokens');
    rightside = regexprep(strtrim(rightsideandcomment{1}{1}),'\(0\)','');
    if length(rightsideandcomment) == 2,
        comment = strtrim(rightsideandcomment{2}{1});
    else
        comment = '';
    end
    % check if state initial condition or if parameter definition
    if isempty(strfind(leftside,'(0)')),
        parametersettings(end+1).name = leftside;
        parametersettings(end).formula = rightside;
        parametersettings(end).notes = comment;
    else
        initialconditions(end+1).name = regexprep(leftside,'\(0\)','');
        initialconditions(end).formula = rightside;
        initialconditions(end).notes = comment;
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [parameterchanges, error] = getParameterChanges(parameterchangesString)
error = '';
% create empty structure
parameterchanges = struct('name',{},'formula',{},'notes',{});
% parse the strings and fill the structures
parameterchangesString = char([double(parameterchangesString) 10]); % need to add a '\n' at the end
allelements = regexp(parameterchangesString,'([^\n]*)=([^\n]*)\n','tokens');
for k = 1:length(allelements),
    element = allelements{k};
    leftside = strtrim(removeCharacters(element{1}));
    rightsideandcomment = strtrim(removeCharacters(element{2}));
    % parse right side and comment
    rightsideandcomment = regexp(rightsideandcomment,'([^%]*)','tokens');
    rightside = strtrim(removeCharacters(rightsideandcomment{1}{1}));
    if length(rightsideandcomment) == 2,
        comment = strtrim(removeCharacters(rightsideandcomment{2}{1}));
    else
        comment = '';
    end
    % check if a piecewisestatement is present
    if ~isempty(strfind(rightside,'{')),
        % construct piecewise statement
        rightside = rightside(2:end-1);  % take away the outer parentheses
        terms = explodePCSB(rightside,',');
        if mod(length(terms),2) == 1, % (if uneven)
            error = sprintf('%sParameter change for ''%s'' does not have a default value defined.',error,leftside);
            return
        end
        if length(terms) <= 2,
            error = sprintf('Wrong definition of parameter setting for parameter ''%s''.\n',leftside);
            return
        end
        formula = 'piecewiseSB(';
        times = terms(1:2:length(terms));
        values = terms(2:2:length(terms));
        triggers = {};
        for k2 = 1:length(times)-1,
            triggers{k2} = sprintf('and(ge(time,%s),le(time,%s))',times{k2},times{k2+1});
        end
        valuestimes = {};
        valuestimes(1:2:length(triggers)+length(values)) = values;
        valuestimes(2:2:length(triggers)+length(values)) = triggers;
        % construct piecewise statement
        formula = 'piecewiseSB(';
        for k2 = 1:length(valuestimes),
            formula = sprintf('%s%s,',formula,valuestimes{k2});
        end
        formula =sprintf('%s)',formula(1:end-1));
    else
        % no piecewise statement detected => just take formula as it is
        formula = rightside;
    end
    % add data
    parameterchanges(end+1).name = leftside;
    parameterchanges(end).formula = formula;
    parameterchanges(end).notes = comment;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [stateevents, error] = getStateEvents(stateeventsString)
error = '';
% create empty structure
eventassignment = struct('variable',{},'formula',{});
stateevents = struct('name',{},'trigger',{},'assignment',eventassignment,'notes',{});
% parse the strings and fill the structures
stateeventsString = char([double(stateeventsString) 13 10]); % need to add a '\n' at the end
if isempty(strtrim(stateeventsString)),
    stateeventsString = strtrim(stateeventsString);
end
allelements = regexp(stateeventsString,'([^\n]*)\n','tokens');
for k=1:length(allelements),
    stateevents(end+1).name = sprintf('StateChange_%d',k);
    % determine the comment
    element = regexp(allelements{k}{1},'([^%]*)','tokens');
    eventdata = element{1}{1};
    if length(element) == 1,
        comment = '';
    else
        comment = element{2}{1};
    end
    stateevents(end).notes = strtrim(removeCharacters(comment));
    % get event data
    terms = explodePCSB(eventdata,',');
    if length(terms) < 2,
        error = sprintf('%s\nState change nr. %d wrongly defined.',error,k);
        return
    end
    % the first term needs to define the time of the event. format: time = xxx
    test = explodePCSB(terms{1},'=');
    if isempty(strcmp(test{1},'time')) || length(test) ~= 2,
        error = sprintf('%s\nState change nr. %d wrongly defined ("time = numeric value" needs to be the first element).',error,k);
        return
    end
    % get trigger function
    stateevents(end).trigger = sprintf('ge(time,%s)',test{2});
    % determine event assignments
    for k2 = 2:length(terms),
        assterms = explodePCSB(terms{k2},'=');
        stateevents(end).assignment(k2-1).variable = assterms{1};
        stateevents(end).assignment(k2-1).formula = assterms{2};
    end
end
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [output] = removeCharacters(input)
% delete all line breaks and tabs from the input string
temp = double(input);
temp(find(temp==13)) = 32;  % replace '\cr' by white space
temp(find(temp==10)) = 32;  % replace '\n' by white space
temp(find(temp==9)) = 32;   % replace '\t' by white space
output = char(temp);
% remove all spaces
%    output = strrep(output,' ','');
return

