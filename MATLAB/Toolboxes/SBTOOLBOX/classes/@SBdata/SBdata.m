function [output] = SBdata(varargin)
% SBdata: creates an object that is intended to contain 
% measurement data - from real experiments or insilico computations.
%
% USAGE:
% ======
% [output] = SBdata()                 creates an empty 
%                                     SBdata object
% [output] = SBdata(struct)           creates an SBdata 
%                                     object from a MATLAB structure 
%                                     in the internal data object format
% [output] = SBdata(inputobject)      construction from a given 
%                                     SBdata object 
% [output] = SBdata('filename.xls')   converting an excel file to 
%                                     an SBdata object
% [output] = SBdata('filename.csv')   converting an CSV (comma
%                                     separated value) file to  
%                                     an SBdata object
% [output] = SBdata(name,notes,time,measurements,componentnames,stimulusFlags) 
%                                     creating an 
%                                     SBdata object 
%                                     from time-series and 
%                                     component names
%
% name: Name for the experiment/measured data
% notes: Notes about the experiment and data
% time: scalar value determining the sampling time for all measurements
%   or time vector of same length as measurement vectors
% measurements: matrix containing the measurements of the components. 
%   One row per timeinstant and one column per component
% componentnames: cell-array with the names of the measured components
% stimulusFlags: vector of same length as componentnames. =0 is indicating
%   a measurement, =1 is indicating a stimulus in the experiment, where the
%   data are collected.
%
% Output Arguments:
% =================
% output: SBdata object 

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
% CHECK THE TYPE OF THE INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==0,
    inputType = 'empty';
elseif nargin == 1,
    if strcmp('SBdata',class(varargin{1})),
        inputType = 'SBdata';
        datainput = varargin{1};
    elseif strcmp('struct',class(varargin{1})),
        inputType = 'datastructure';
        datastructure = varargin{1};
    elseif strcmp('char',class(varargin{1})),
        % assume filenamec given as input argument
        % check if '.xls' given as extension. If yes, then import data from
        % excel file format
        filename = varargin{1};
        if ~isempty(strfind(filename,'.xls')),
            inputType = 'XLSdataFile';
        elseif ~isempty(strfind(filename,'.csv')),
            inputType = 'CSVdataFile';
        else
            error('Unknown filetype!');
        end
    else 
        error('Input argument of unknown type');
    end
elseif nargin == 6,
    inputType = 'simple';
    name = varargin{1};
    notes = varargin{2};
    time = varargin{3};
    measurements = varargin{4};
    componentnames = varargin{5};
    stimulusFlags = varargin{6};
    % check the input for consistency
    if size(measurements,2) ~= length(componentnames), error('Different number of measured components and component names.'); end
    if size(measurements,2) ~= length(stimulusFlags), error('Different number of component names and stimuli flags.'); end
else
    error('Wrong number of input arguments.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCT THE SBMODEL OBJECT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp('empty',inputType),
    % Create empty data structure
    % measurements substructure
    measurementsStruct = struct('name',{},'stimulusflag',{},'data',{},'noiseoffset',{},'noisevariance',{});
    % Create data structure
    datastructure = struct('name','','notes','','time',{},'measurements',measurementsStruct);
    % construct the model object
    datastructure(1).name = 'untitled';
    datastructure(1).notes = 'no notes';
    datastructure(1).time = [];
    datastructure(1).measurements = measurementsStruct;
    output = class(datastructure,'SBdata');
elseif strcmp('SBdata',inputType),
    % copy the model object
    output = datainput;
elseif strcmp('datastructure',inputType),
    % check if the given structure is a data structure (only check the
    % top-level fields)
    checkFields = {'name','notes','time','measurements'};
    for k = 1:length(checkFields),
        if ~isfield(datastructure,checkFields{k}),
            errorMsg = sprintf('Given structure is not a valid internal SBdata structure.');
            error(errorMsg);
        end
    end
    % construct the data object
    output = class(datastructure,'SBdata');
elseif strcmp('XLSdataFile',inputType),
    % check if a file with given filename exists
    [path,filename,ext,versn] = fileparts(filename);
    filename = strcat(filename,'.xls');
    if ~exist(filename),
        errorMsg = sprintf('XLS data file does not exist in current directory.');
        error(errorMsg);
    end
    % If file exists then import it
    [datastructure,errorMsg] = SBimportXLSdata(filename);
    % Check if error occurred while importing the XLS data
    if length(errorMsg) ~= 0,
        error(errorMsg);
    end
    % construct the data object
    output = class(datastructure,'SBdata');
elseif strcmp('CSVdataFile',inputType),
    % check if a file with given filename exists
    [path,filename,ext,versn] = fileparts(filename);
    filename = strcat(filename,'.csv');
    if ~exist(filename),
        errorMsg = sprintf('CSV data file does not exist in current directory.');
        error(errorMsg);
    end
    % If file exists then import it
    [datastructure,errorMsg] = SBimportCSVdata(filename);
    % Check if error occurred while importing the CSV data
    if length(errorMsg) ~= 0,
        error(errorMsg);
    end
    % construct the data object
    output = class(datastructure,'SBdata');
elseif strcmp('simple',inputType);
    % create datastructure from simple input arguments
    datastructure.name = name;
    datastructure.notes = notes;
    % check if time given as deltaT or time vector
    if length(time) == 1,
        time = 0:time:time*(size(measurements,1)-1);
    end
    % check length of time vector
    if length(time) ~= size(measurements,1),
        error('Length of given time vector does not match length of measurements.');
    end
    % check if time vector is monotonically increasing
    test = time(2:end)-time(1:end-1);
    if ~isempty(find(test <= 0))
        error('Given time vector is not monotonically increasing.');
    end
    datastructure.time = time;
    % cycle through all measured components and add info about them
    for k = 1:length(componentnames),
        datastructure.measurements(k).name = componentnames{k};
        datastructure.measurements(k).stimulusflag = stimulusFlags(k);
        datastructure.measurements(k).data = measurements(:,k);
        datastructure.measurements(k).noiseoffset = 0;
        datastructure.measurements(k).noisevariance = 0;
    end
    % construct the data object
    output = class(datastructure,'SBdata');
end
return
