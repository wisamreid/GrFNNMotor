function [modelFunctionNames] = getSBmodelFunctionNames(model)
% getSBmodelFunctionNames
% gives back a structure containing all function names defined
% in this SBmodel
%
% USAGE:
% ======
% [modelFunctionNames] = getSBmodelFunctionNames(SBmodel) 
%
% SBmodel: SBmodel 
%
% modelFunctionNames: structure with all function names

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

% CHECK IF SBmodel
if ~strcmp('SBmodel',class(model)),
    error('Function only defined for SBmodels.');
end
% get the datastructure of the model
sbm = SBstruct(model);

modelFunctionNames = {};
% fetch all names of defined function within given SBmodel
for index = 1 : length(sbm.functions),
    modelFunctionNames{index} = sbm.functions(index).name;
end

return