function [output] = SBinitialconditions(varargin)
% SBinitialconditions: Returns a column vector with the values of the 
% states in the SBmodel as elements. (The values are the initial conditions
% stored in the model). The function can also be used to set initial
% conditions of the model, but only if the model is given as an SBmodel and
% not as an m-file.
%
% USAGE:
% ======
% [model] = SBinitialconditions(model,values) 
% [values] = SBinitialconditions(model,statenames) 
% [values] = SBinitialconditions(model)
%
% model: SBmodel, or ODE file model description (in the latter case the
%   initial conditions can't be set!)
% statenames: cell-array with statenames for which to give back the initial
%   conditions.
% values: vector with initial conditions
%
% Output Arguments:
% =================
% model: SBmodel with changed initial conditions

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

if nargin == 1,
    % one input argument => assume it is a model description 
    % (SBmodel or ODE file name)
    model = varargin{1};
    % obtain the initial conditions
    if strcmp('SBmodel',class(model)),
        SBstructure = SBstruct(model);
        values = [SBstructure.states.initialCondition];
    elseif strcmp('char',class(model)),
        values = feval(model);
    else
        error('Unknown input type');
    end
    % return the values
    output = values(:);
elseif nargin == 2,    
    % two input arguments => assume the first is a model description 
    model = varargin{1};
    % check second input argument
    if strcmp(class(varargin{2}),'double'),
        % the second argument is a double ... set new initial conditions.
        if ~strcmp('SBmodel',class(model)),
            error(sprintf('Changing initial conditions works only for SBmodels, not for ODE file\ndescriptions of models.'));
        end
        SBstructure = SBstruct(model);
        values = varargin{2};
        if length(values) ~= length(SBstructure.states),
            error('The length of the values vector does not match the number of states in the model.');
        end
        % Set new initial values
        cellvalues = mat2cell(values(:)',1,ones(1,length(values)));
        [SBstructure.states.initialCondition] = deal(cellvalues{:});
        % Return the SBmodel with new initial values
        output = SBmodel(SBstructure);
    else
        % second input argument are statenames for which to return initial
        % conditions
        statenames = varargin{2};
        [allnames, dummy, initialconditions] = SBstates(model);
        if strcmp(class(statenames),'char'), 
            statenames = {statenames}; 
        end
        output = [];
        for k=1:length(statenames),
            index = strmatch(statenames{k},allnames,'exact');
            if isempty(index),
                error(sprintf('State ''%s'' does not exist in the model.',statenames{k}));
            end
            output(k) = initialconditions(index);
        end
    end
    
else
    error('Incorrect number of input arguments.');
end

