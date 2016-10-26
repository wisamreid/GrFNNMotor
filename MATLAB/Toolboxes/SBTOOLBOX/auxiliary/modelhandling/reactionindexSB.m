function [output] = reactionindexSB(model,reactionname)
% reactionindexSB: returns the number of the reaction 'reactionname' in model
% 'model'. If the reaction does not exist then [] is returned.
%
% Output Arguments:
% =================
% output = index of the reaction 'reactionname' in the model.
%          If 'reactionname' is not a reaction in the model, [] is returned.

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

output = strmatch(reactionname,SBreactions(model),'exact');
if isempty(output),
    output = [];
end
return