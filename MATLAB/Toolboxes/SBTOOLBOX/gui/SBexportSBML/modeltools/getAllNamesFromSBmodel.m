function [nameList] = getAllNamesFromSBmodel(model, varargin)
% getAllNamesFromSBmodel
% collects the names of all components used in the given SBmodel
%
%
% USAGE:
% ======
% [nameList] = getAllNamesFromSBmodel(model)
%
% model: SBmodel
% 
% nameList: list of component names used in the model
%

% Information:
% ============
% Systems Biology Toolbox for MATLAB
% Copyright (C) 2005-2007 Fraunhofer-Chalmers Research Centre
% Author of the toolbox: Henning Schmidt, henning@hschmidt.de
% Programmed by Gunnar Drews, Student at University of 
% Rostock Department of Computerscience, gunnar.drews@uni-rostock.de
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
% CHECK IF SBmodel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp('SBmodel',class(model)),
    error('Function only defined for SBmodels.');
end

silentFlag = 0;
if nargin == 1,
    silentFlag = 0;
elseif nargin == 2,
    silentFlag = varargin{1};
end

componentNames = [];
nameIndex = 1;

% get model structure and count model components
sbm = SBstruct(model);
statesCount = length(sbm.states);
parametersCount = length(sbm.parameters);
variablesCount = length(sbm.variables);
reactionsCount = length(sbm.reactions);
eventsCount = length(sbm.events);
functionsCount = length(sbm.functions);

% fetch state names
for k1 = 1 : statesCount,
    componentNames{nameIndex} = sbm.states(k1).name;
    nameIndex = nameIndex + 1;
end

% fetch parameter names
for k1 = 1 : parametersCount,
    componentNames{nameIndex} = sbm.parameters(k1).name;
    nameIndex = nameIndex + 1;
end

% fetch variable names
for k1 = 1 : variablesCount,
    componentNames{nameIndex} = sbm.variables(k1).name;
    nameIndex = nameIndex + 1;
end

% fetch reaction names
for k1 = 1 : reactionsCount,
    componentNames{nameIndex} = sbm.reactions(k1).name;
    nameIndex = nameIndex + 1;
end

% fetch event names
for k1 = 1 : eventsCount,
    componentNames{nameIndex} = sbm.events(k1).name;
    nameIndex = nameIndex + 1;
end

% fetch functions names
for k1 = 1 : functionsCount,
    componentNames{nameIndex} = sbm.functions(k1).name;
    nameIndex = nameIndex + 1;
end

nameList = componentNames;
return