function [varargout] = SBsimulate(varargin)
% SBsimulate: simulates an SBmodel or an ODE file created by
% SBcreateODEfile. In case the model is given as an SBmodel, not only the 
% trajectories for the models states are determined, but also the time 
% trajectories for all variables and reaction rates in the model 
% 
% The SBsimulate function can deal automatically with models that contain 
% state events. However, this is only possible in the case where a model is
% specified as an SBmodel, not as an ODE file model.
%
% USAGE:
% ======
% [output] = SBsimulate(model)         
% [output] = SBsimulate(model,time)         
% [output] = SBsimulate(model,method,time)         
% [output] = SBsimulate(model,method,time,ic)         
% [output] = SBsimulate(model,method,time,ic,options)         
%
% model: SBmodel or ODE file model description
% time: scalar value for simulation end time (default start time is 0),
%       vector with two elements indicating start and end time, or vector
%       with time instants for which the simulation results are to be
%       returned.
% method: name of a MATLAB integrator (ode45, ode23s, etc.) 
% ic: initial conditions of the model to be simulated
% options: standard MATLAB integrator options, set by ODESET
%
% DEFAULT VALUES:
% ===============
% time: 20 time units (=> tspan: [0 20])
% method: 'ode23s'
% ic: the initial conditions stored in the model
% options: []
%
% Output Arguments:
% =================
% If no output arguments are given, the result of the simulation is plotted
% 1) online during the simulation for monitoring purposes - the simulation 
% can here be stopped at every instant by just clicking on the "Stop"
% button. 2) after the simulation is finished the simulation data is plotted 
% using the SBplot function, allowing to browse the data.
%
% The output of SBsimulate is realized as a structure:
% output.time:              vector with time instants for the results
% output.states:            cell-array with state names
% output.statevalues:       matrix with state values. Each row corresponds to
%                           one time instant.
%
% The following fields are only present in case the model as been given
% as an SBmodel:
%
% output.variables:         cell-array with variable names
% output.variablevalues:    matrix with variable values. Each row corresponds to
%                           one time instant.
% output.reactions:         cell-array with reaction names
% output.reactionvalues:    matrix with reaction values. Each row corresponds to
%                           one time instant.

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
% INITIALIZE STUFF
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = [];
x = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DEAL WITH VARIABLE ARGUMENT NUMBER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 1,
    model = varargin{1};
    method = 'ode23s';
    tspan = [0 20];
    ic = [];
    options = [];
elseif nargin == 2,
    model = varargin{1};
    method = 'ode23s';
    if length(varargin{2}) == 1,
        tspan = [0 varargin{2}];
    else
        tspan = varargin{2};
    end
    ic = [];
    options = [];
elseif nargin == 3,
    model = varargin{1};
    method = varargin{2};
    if length(varargin{3}) == 1,
        tspan = [0 varargin{3}];
    else
        tspan = varargin{3};
    end
    ic = [];
    options = [];
elseif nargin == 4,
    model = varargin{1};
    method = varargin{2};
    if length(varargin{3}) == 1,
        tspan = [0 varargin{3}];
    else
        tspan = varargin{3};
    end
    ic = varargin{4};
    options = [];
elseif nargin == 5,
    model = varargin{1};
    method = varargin{2};
    if length(varargin{3}) == 1,
        tspan = [0 varargin{3}];
    else
        tspan = varargin{3};
    end
    ic = varargin{4};
    options = varargin{5};
else
    error('Incorrect number of input arguments');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK NUMBER OF OUTPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout ~= 0 && nargout ~= 1,
    error('Incorrect number of output arguments');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS SBMODEL OR ODE FILE - CHECK IF EVENTS PRESENT IN SBmodel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
eventFlag = 0;
if strcmp('SBmodel',class(model)),
    % check for events
    if ~isempty(SBevents(model)) > 0,
        eventFlag = 1;
    end
    [ODEfctname, ODEfilefullpath, DATAfctname, EVENTfctname, EVENTASSGNfctname] = SBcreateTempODEfile(model,1,eventFlag);
else
    % ODEfctname of ODE file
    ODEfctname = model;
    % Check if file exists
    if exist(strcat(ODEfctname,'.m')) ~= 2,
        error('ODE file could not be found.');
    end
end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% UPDATE INITIAL CONDITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(ic),
    ic = feval(ODEfctname);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% IF NO OUTPUT ARGUMENTS GIVEN THEN DISPLAY ONLINE INFORMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout == 0,
    figure; drawnow;
    options = odeset(odeset(options),'OutputFcn',@odeplot);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SIMULATE MODEL WITH GIVEN METHOD (DIFFERENT HANDLING IN CASE OF EVENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if eventFlag == 0,
    [t,x] = feval(method,ODEfctname,tspan,ic,options);
else
    % handle the events that are present in the model
    [t,x] = eventsimulateSB(ODEfctname,method,tspan,ic,options,EVENTfctname,EVENTASSGNfctname);
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLOSE ONLINE PLOT AGAIN IF IT WAS PRESENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout == 0,
    close gcf;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DETERMINE OUTPUT VARIABLE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
output.time = t;
output.states = SBstates(ODEfctname)';
output.statevalues = x;
if strcmp('SBmodel',class(model)),
    % determine the simulation data using the 'DATAfctname'
    [dataoutput] = feval(DATAfctname, t,x);
    output.variables = dataoutput.variables;
    output.variablevalues = dataoutput.variablevalues;
    output.reactions = dataoutput.reactions;
    output.reactionvalues = dataoutput.reactionvalues;
end
if nargout == 0,
    % prepare data for plotting
    time = output.time;
    datanames = {};
    dataindex = 1;
    for k = 1:length(output.states),
        datanames{dataindex} = sprintf('%s (state)',output.states{k});
        dataindex = dataindex + 1;
    end
    if strcmp('SBmodel',class(model)),
    for k = 1:length(output.variables),
        datanames{dataindex} = sprintf('%s (variable)',output.variables{k});
        dataindex = dataindex + 1;
    end
    for k = 1:length(output.reactions),
        datanames{dataindex} = sprintf('%s (reaction rate)',output.reactions{k});
        dataindex = dataindex + 1;
    end
    datavalues = [output.statevalues, output.variablevalues, output.reactionvalues];
    else
        datavalues = [output.statevalues];
    end
    SBplot(time,datavalues,datanames);
elseif nargout == 1,
    varargout{1} = output;
else 
    error('Incorrect number of output arguments!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DELETE FILE IF SBMODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp('SBmodel',class(model)),
    deleteTempODEfileSB(ODEfilefullpath);
end

