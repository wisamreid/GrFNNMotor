function [output] = isstateSB(model,name)
% isstateSB: checks if "name" is a state in the provided model.
% This function works both for SBmodels and ODE files. The check is of
% course case sensitive
%
% Output Arguments:
% =================
% output: =1 if "name" is a state, =0 if not

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

% get all states of the model
allStates = SBstates(model);

% check if "name" is a state.
output = 0;
for k = 1:length(allStates),
    if strcmp(strtrim(name),allStates{k}),
        output = 1;
        break;
    end
end
return