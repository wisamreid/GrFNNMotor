function [boolvalue] = isSBmodelName(name, varargin)
% isSBmodelName
% checks wether a given name complies with the used name standard
% of the SBmodel
%
%
% USAGE:
% ======
% [boolvalue] = isSBmodelName(name)
%
% name: string to test SBmodel compliance
% 
% boolvalue: true if the given name complies with the SBmodel name
%            conventions 
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

silentFlag = 0;
if (nargin == 2)
    silentFlag = 0;
elseif (nargin == 3)
    silentFlag = varargin{1};
end

boolvalue = false;

%use regular expression to test, wether given name complies with SBmodel
%name conventions
result=regexp(name, '([A-Z]|[a-z]|_)+([A-Z]|[a-z]|[0-9]|_)*');
if (result == 1)
    boolvalue = true;
end

return