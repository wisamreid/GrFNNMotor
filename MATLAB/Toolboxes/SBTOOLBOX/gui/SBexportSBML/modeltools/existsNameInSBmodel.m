function [boolvalue] = existsNameInSBmodel(model, name, varargin)
% existsNameInSBmodel
% checks wether a given name is already used in the given SBmodel
%
%
% USAGE:
% ======
% [boolvalue] = existsNameInSBmodel(model, name)
%
% model: SBmodel
% name: string to test model for existence
% 
% boolvalue: true if there is already defined a component with the given
%            name otherwise false
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
if (nargin == 2)
    silentFlag = 0;
elseif (nargin == 3)
    silentFlag = varargin{1};
end

boolvalue = false;

% fetch all names used within given SBmodel and test wether given name is
% already present
names = getAllNamesFromSBmodel(model);
rowPos = matchStringOnArray(name, names);
if (rowPos > 0)
    boolvalue = true;
end

return