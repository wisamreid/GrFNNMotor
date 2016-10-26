function [expmodel] = SBmergemodexp(model, experiment)
% SBmergemodexp: combines an experiment with a model, and produces a "merged 
% model", as output. The original model and experiment arguments can be
% given as textfiles, or SBexp objects. The output model is an
% SBmodel.
%
% DESCRIPTIONS + SYNTAX:
% ======================
% To fill in!
%
% USAGE:
% ======
% [expmodel] = SBmergemodexp(model, experiment)        
% [expmodel] = SBmergemodexp(model, experimentfile)        
%
% model: SBmodel 
% experiment: An SBexp object describing an experiment that should be done
%             with the model
% experimentfile: String with the name of an experiment file
%
% Output Arguments:
% =================
% Merged model containing the original model and the experimental settings.

% Information:
% ============
% Systems Biology Toolbox for MATLAB
% Copyright (C) 2005-2007 Fraunhofer-Chalmers Research Centre
% Author of the toolbox: Henning Schmidt, henning@hschmidt.de
% Functions authors: Joachim Almquist and Nils Lämås, FCC 
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Handle input arguments 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp(class(model),'SBmodel'),
    error('The first input argument needs to be an SBmodel.');
else
    m = SBstruct(model);
end
if strcmp(class(experiment),'SBexp'),
    e = SBstruct(experiment);
elseif strcmp(class(experiment),'char'),
    if ~exist(experiment),
        error('Experiment file does not exist in current directory.');
    else
        e = SBstruct(SBexp(experiment));
    end
else
    error('Incorrect second (experiment) input argument.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize the output model structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
o = m;                    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% find all parameters and initial conditions listed in the experiment 
% and erase them in the output model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(e.parametersettings)
    try
        index = find_index({o.parameters.name}, e.parametersettings(i).name);
        o.parameters(index).value = [];
    catch
        o.parameters(end+1).name = e.parametersettings(i).name;
    end
end

for i = 1:length(e.initialconditions)
    try
        index = find_index({o.states.name}, e.initialconditions(i).name);
        o.states(index).initialCondition = [];
    catch
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% evaluate the remaining parameters and initial conditions in the original model 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(o.parameters)
    if ~isempty(o.parameters(i).value)
        eval([o.parameters(i).name '=' num2str(o.parameters(i).value) ';']);
    end
end

for i = 1:length(o.states)
    if ~isempty(o.states(i).initialCondition)
        eval([o.states(i).name '=' num2str(o.states(i).initialCondition) ';']);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define parametersettings and initial conditions in the output model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
count = 1;
while count > 0
    count = 0;
    for i = 1:length(e.parametersettings)
        index = find_index({o.parameters.name}, e.parametersettings(i).name);
        if ~isempty(e.parametersettings(i).formula) && isempty(o.parameters(index).value)
            try
                count = count+1;
                o.parameters(index).value = eval(e.parametersettings(i).formula);
                eval([o.parameters(index).name '=' num2str(o.parameters(index).value) ';']);              
            catch
                count = count-1;
            end
        end
    end
    
    for i = 1:length(e.initialconditions)
        index = find_index({o.states.name}, e.initialconditions(i).name);
        if ~isempty(e.initialconditions(i).formula) && isempty(o.states(index).initialCondition)
            try
                count = count+1;
                o.states(index).initialCondition = eval(e.initialconditions(i).formula);
                eval([o.states(index).name '=' num2str(o.states(index).initialCondition) ';']);               
            catch                
                count = count-1;
            end
        end
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define parameter changes in the output model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save old variables
oldvariables = o.variables;
o.variables = [];
% add new variables as first in the list
for i = 1:length(e.parameterchanges)    
    try        
        index = find_index({e.parametersettings.name}, e.parameterchanges(i).name);
        disp('Error. At least one parameter is defined in parameter settings AND in parameter changes.')
        return
    catch
    end
    try
        index = find_index({o.parameters.name}, e.parameterchanges(i).name);
        o.parameters(index) = [];
    catch
    end
    try
        index = find_index({o.variables.name}, e.parameterchanges(i).name);
        o.variables(index) = [];
    catch
    end
    o.variables(end+1).name = e.parameterchanges(i).name;
    o.variables(end).formula = e.parameterchanges(i).formula;
    o.variables(end).notes = e.parameterchanges(i).notes;
    o.variables(end).type = '';
    o.variables(end).compartment = '';
    o.variables(end).unittype = '';
end
% add old variables
for k = 1:length(oldvariables),
    o.variables(end+1).name = oldvariables(k).name;
    o.variables(end).formula = oldvariables(k).formula;
    o.variables(end).notes = oldvariables(k).notes;
    o.variables(end).type = oldvariables(k).type;
    o.variables(end).compartment = oldvariables(k).compartment;
    o.variables(end).unittype = oldvariables(k).unittype;    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% define events in the output model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:length(e.stateevents)
    o.events(end+1) = e.stateevents(i);
end
if (length([o.parameters.value]) ~= length(o.parameters)) || (length([o.states.initialCondition]) ~= length(o.states))
    disp('Error. At least one parameter/initial condition has no defined value. Check experiment for inconsistency.')
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% return the experiment model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
expmodel = SBmodel(o);
return


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The "find index" function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [index] = find_index(name_array, comp_string)
i = 0;
l = length(name_array);
while i < l
    i = i+1;
    if strcmp(comp_string, name_array(i))
        index = i;
    end
end
return