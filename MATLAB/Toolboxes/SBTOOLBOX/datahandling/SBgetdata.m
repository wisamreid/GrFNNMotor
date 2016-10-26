function [varargout] = SBgetdata(varargin)
% SBgetdata
% This functions allows to extract information about a data structure
% 
% USAGE:
% ======
% [] = SBgetdata(data)
% [time, componentNamesMeasured, dataMeasured, ...
%       componentNamesStimuli, dataStimuli] = SBgetdata(data)
% [] = SBgetdata(data,componentname)
% [time,datavalues,stimulusflag,noiseoffset,noisevariance,deltaT] = ...
%       SBgetdata(data,componentname)
% [time,datavalues,stimulusflag,noiseoffset,noisevariance,deltaT] = ...
%       SBgetdata(data,componentname,tol)
% 
% data: SBdata object containing the data
% componentname: name of the component that is to be extracted
% tol: tolerance for the decision if equally spaced time steps in time
%      vector (used for the determination of the sampling time). Value given 
%      in percent of the mean of the time steps for the different components.
%
% DEFAULT VALUES:
% ===============
% tol: 1%
%
% Output Arguments:
% =================
% If no output argument is given, the function outputs the information to
% the MATLAB window.
%
% time: time vector of measurement instants for the overall 
%   measurement data and for single components
% componentNamesMeasured: cell-array containing the names of the measured
%   components
% dataMeasured: matrix containing all the measurements of the SBdata
%   object. One row per time instant and one column per measured component.
%   The ordering of the columns corresponds to the ordering of the names in
%   the "componentNamesMeasured" output variable.
% componentNamesStimuli: cell-array containing the names of the stimulated
%   components
% dataStimuli: matrix containing all the stimuli data of the SBdata
%   object. One row per time instant and one column per stimulated component.
%   The ordering of the columns corresponds to the ordering of the names in
%   the "componentNamesStimuli" output variable.
%
% datavalues: vector of measurement data for given component
% stimulusflag: flag indicating if component is a stimulus (=1) 
%   or a measurement (=0)
% noiseoffset: noise offset defined for given component
% noisevariance: noise variance for given component
% deltaT: sampling time for given component. If non equidistantly sampled,
%   0 is returned.
%
% Please note that 0 is default value for noise information in case no 
% noise information is specified for a component.

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
    componentname = 'nonenone';
elseif nargin == 2,
    data = SBstruct(varargin{1});
    componentname = varargin{2};
    tol = 1;
elseif nargin == 3,
    data = SBstruct(varargin{1});
    componentname = varargin{2};
    tol = varargin{3};
else 
    error('Incorrect number of input arguments.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK IF MODEL CONTAINS MEASUREMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(data.measurements) == 0,
    error('The model does not contain any measurements');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET INFORMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp(componentname,'nonenone'),
    % get information about single component
    % first get its index
    foundComponent = 0;
    for k = 1:length(data.measurements),
        if strcmp(componentname,data.measurements(k).name),
            foundComponent = 1;
            break;
        end
    end
    % check if component found
    if foundComponent == 0,
        text = sprintf('Given component not available in the ''%s'' data',data.name);
        error(text);
    end
    % extract stimulus flag
    stimulusflag = data.measurements(k).stimulusflag;
    % extract information about the given component and time info
    datavalues = data.measurements(k).data;
    time = data.time;
    % process NaN values in the data vector
    time = time(find(~isnan(datavalues)));
    datavalues = datavalues(find(~isnan(datavalues)));
    noiseoffset = data.measurements(k).noiseoffset;
    noisevariance = data.measurements(k).noisevariance;
    % detect the sampling time deltaT
    % analyzing the time vector
    deltaTime = time(2:end)-time(1:end-1);
    % determine deltaT (function returns 0 if variation larger than tolerance)
    deltaT = checksamplingtimevectorSB(deltaTime,tol);
else
    % get time vector
    time = data.time;
    % get component and stimuli names
    componentNamesMeasured = {};
    componentNamesStimuli = {};
    dataMeasured = [];
    dataStimuli = [];
    measuredNamesIndex = 1;
    stimuliNamesIndex = 1;
    for k = 1:length(data.measurements),
        if data.measurements(k).stimulusflag == 0,
            componentNamesMeasured{measuredNamesIndex} = data.measurements(k).name;
            dataMeasured(:,measuredNamesIndex) = data.measurements(k).data(:);
            measuredNamesIndex = measuredNamesIndex + 1;
        else
            componentNamesStimuli{stimuliNamesIndex} = data.measurements(k).name;
            dataStimuli(:,stimuliNamesIndex) = data.measurements(k).data(:);
            stimuliNamesIndex = stimuliNamesIndex + 1;
        end
    end
    componentNamesMeasured = componentNamesMeasured(:);
    componentNamesStimuli = componentNamesStimuli(:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS VARIABLE OUTPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout == 0,
    if strcmp(componentname,'nonenone'),
        if length(componentNamesMeasured) > 0,
            text = sprintf('Measured components available in ''%s'' are:\n',data.name);
            for k = 1:length(componentNamesMeasured),
                text = sprintf('%s%s\n',text,componentNamesMeasured{k});
            end 
            disp(text);
        else
            disp('No measured components present');
        end
        if length(componentNamesStimuli) > 0,
            text = sprintf('Stimulus components available in ''%s'' are:\n',data.name);
            for k = 1:length(componentNamesStimuli),
                text = sprintf('%s%s\n',text,componentNamesStimuli{k});
            end 
            disp(text);
        else
            disp('No stimulus components present');
        end
    else
        text = sprintf('Information about component ''%s'' in ''%s'':\n',componentname,data.name);
        if stimulusflag == 0,
            text = sprintf('%sMeasured component\n',text);
        else
            text = sprintf('%sStimulus component\n',text);
        end
        text = sprintf('%sStart time: %g\n',text,min(time));
        text = sprintf('%sEnd time: %g\n',text,max(time));
        if deltaT == 0,
            text = sprintf('%sNo equidistant sampling\n',text);
        else
            text = sprintf('%sSampling time: %g\n',text,deltaT);
        end
        text = sprintf('%sNumber of samples: %d\n',text,length(datavalues));
        text = sprintf('%sNoise offset: %g\n',text,noiseoffset);
        text = sprintf('%sNoise variance: %g\n',text,noisevariance);
        disp(text);
    end
elseif nargout == 1,
    varargout{1} = time;
elseif nargout == 2,
    varargout{1} = time;
    if strcmp(componentname,'nonenone'),
        varargout{2} = componentNamesMeasured;
    else
        varargout{2} = datavalues;
    end
elseif nargout == 3,
    varargout{1} = time;
    if strcmp(componentname,'nonenone'),
        varargout{2} = componentNamesMeasured;
        varargout{3} = dataMeasured;
    else
        varargout{2} = datavalues;
        varargout{3} = stimulusflag;
    end
elseif nargout == 4,
    varargout{1} = time;
    if strcmp(componentname,'nonenone'),
        varargout{2} = componentNamesMeasured;
        varargout{3} = dataMeasured;
        varargout{4} = componentNamesStimuli;
    else
        varargout{2} = datavalues;
        varargout{3} = stimulusflag;
        varargout{4} = noiseoffset;
    end
elseif nargout == 5,
    varargout{1} = time;
    if strcmp(componentname,'nonenone'),
        varargout{2} = componentNamesMeasured;
        varargout{3} = dataMeasured;
        varargout{4} = componentNamesStimuli;
        varargout{5} = dataStimuli;
    else
        varargout{2} = datavalues;
        varargout{3} = stimulusflag;
        varargout{4} = noiseoffset;
        varargout{5} = noisevariance;
    end
elseif nargout == 6,
    varargout{1} = time;
    if strcmp(componentname,'nonenone'),
        error('Incorrect number of output arguments.');
    else
        varargout{2} = datavalues;
        varargout{3} = stimulusflag;
        varargout{4} = noiseoffset;
        varargout{5} = noisevariance;
        varargout{6} = deltaT;
    end
else
    error('Incorrect number of output arguments.');
end
return

