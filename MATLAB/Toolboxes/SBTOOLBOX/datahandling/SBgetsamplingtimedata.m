 function [deltaT, componentNamesMeasured, componentNamesStimuli, ...
     componentDeltaTmeasured, componentDeltaTstimuli] = SBgetsamplingtimedata(varargin)
% SBgetsamplingtimedata
% This function determines the sampling times that are used in a given
% SBdata object. This is done by checking if the sampling instants are
% equally spaced for each measured component. 
% 
% USAGE:
% ======
% [deltaT, componentNamesMeasured, componentNamesStimuli, ...
%    componentDeltaTmeasured, componentDeltaTstimuli] = SBgetsamplingtimedata(data)
% [deltaT, componentNamesMeasured, componentNamesStimuli, ...
%    componentDeltaTmeasured, componentDeltaTstimuli] = SBgetsamplingtimedata(data,tol)
%
% data: SBdata object containing the data
% tol: tolerance for the decision if equally spaced time steps in time
%      vector. Value given in percent of the mean of the time steps for the
%      different components.
%
% DEFAULT VALUES:
% ===============
% tol: 1%
%
% Output Arguments:
% =================
% deltaT: If the same sampling time is used for all measured components it 
%   is returned in this variable. If different sampling times are used or if
%   non equidistant sampling is present, 0 is returned.
% componentNamesMeasured: cell-array containing the names of the measured
%   components present in the data object.
% componentNamesStimuli: cell-array containing the names of the stimulated
%   components present in the data object.
% componentDeltaTmeasured: vector of sampling times for the different
%   measured components (same order as componentNamesMeasured). 0 indicates non
%   equidistant sampling. 
% componentDeltaTstimuli: vector of sampling times for the different
%   stimulated components (same order as componentNamesStimuli). 0 indicates non
%   equidistant sampling. 

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
% PROCESS VARIABLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 1,
    data = SBstruct(varargin{1});
    tol = 1;
elseif nargin == 2,
    data = SBstruct(varargin{1});
    tol = varargin{2};
else
    error('Incorrect numer of input arguments.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK IF MODEL CONTAINS MEASUREMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(data.measurements) == 0,
    error('The model does not contain any measurements');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize output data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deltaT = [];
componentDeltaTmeasured = [];
componentDeltaTstimuli = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get all componentnames
[time,componentNamesMeasured,dataMeasured,componentNamesStimuli,dataStimuli] = SBgetdata(SBdata(data));
% cycle through componentnames and get sampling time 
for k = 1:length(componentNamesMeasured),
    [a,b,b,d,e,deltaTk] = SBgetdata(SBdata(data),componentNamesMeasured{k},tol);
    componentDeltaTmeasured = [componentDeltaTmeasured, deltaTk];
end
for k = 1:length(componentNamesStimuli),
    [a,b,b,d,e,deltaTk] = SBgetdata(SBdata(data),componentNamesStimuli{k},tol);
    componentDeltaTstimuli = [componentDeltaTstimuli, deltaTk];
end
% determine deltaT (function returns 0 if variation larger than tolerance)
deltaT = checksamplingtimevectorSB([componentDeltaTmeasured componentDeltaTstimuli],tol);
return




