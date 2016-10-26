function [names,formulas,reversibleFlag,varargout] = SBreactions(model,varargin)
% SBreactions: Returns information about the reactions in an SBmodel.
% If a state vector is passed additionally, the corresponding reaction
% rates are returned also.
%
% USAGE:
% ======
% [names,formulas,reversibleFlag] = SBreactions(model)
% [names,formulas,reversibleFlag,reactionrates] = SBreactions(model,statevector)
% [names,formulas,reversibleFlag,reactionrates] = SBreactions(model,time,statevector)
%
% model: SBmodel (function can not be used on M-file model)
% statevector: vector with corresponding statevalues
% time: time instant of the statevector (only needed for time variant
% systems)
%
% DEFAULT VALUES:
% ===============
% statevector: not needed 
% time: 0  (makes no difference for time invariant systems)
%
% Output Arguments:
% =================
% names: cell-array with models reaction names
% formulas: cell-array with formuas for the kinetic laws of the reactions
% reversibleFlag: vector with same number of entries as reactions. An entry
%   of 0 indicates an irreversible reaction, an entry of 1 indicates a
%   reversible reaction. The ordering is the same as the reaction names.
% reactionrates: the reaction rates at the given state and time

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
    names = {sbm.reactions.name};
    formulas = {sbm.reactions.formula};
    reversibleFlag = [sbm.reactions.reversible];
else
    error('The function can only be used on SBmodels, not on M-file ODE models');
end
names = names(:);
formulas = formulas(:);
reversibleFlag = reversibleFlag(:);

% check if the reaction rates need to be determined:
if nargin == 2,    
    statevector = varargin{1};
    time = 0;
    % create data file (using SBcreateTempODEfile function)
    [ODEfctname, ODEfilefullpath, DATAfctname] = SBcreateTempODEfile(model,1);
    data = feval(DATAfctname, time, statevector);
    % delete all temporary m files
    deleteTempODEfileSB(ODEfilefullpath);
    varargout{1} = data.reactionvalues';
elseif nargin == 3,
    time = varargin{1};
    statevector = varargin{2};
    % create data file (using SBcreateTempODEfile function)
    [ODEfctname, ODEfilefullpath, DATAfctname] = SBcreateTempODEfile(model,1);
    data = feval(DATAfctname, time, statevector);
    % delete all temporary m files
    deleteTempODEfileSB(ODEfilefullpath);
    varargout{1} = data.reactionvalues';
end
return