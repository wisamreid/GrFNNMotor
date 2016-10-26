function [functionsMATLAB] = SBfunctionsMATLAB(model)
% SBfunctionsMATLAB: Returns information about the MATLAB functions in an
% SBmodel. 
%
% USAGE:
% ======
% [functionsMATLAB] = SBfunctionsMATLAB(model)
%
% model: SBmodel (function can not be used on M-file model)
%
% Output Arguments:
% =================
% functionsMATLAB: string containing eventual MATLAB functions defined in
%   the model.

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS SBMODEL OR ODE FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp('SBmodel',class(model)),
    sbm = SBstruct(model);
    functionsMATLAB = sbm.functionsMATLAB;
else
    error('The function can only be used on SBmodels, not on M-file ODE models');
end
return