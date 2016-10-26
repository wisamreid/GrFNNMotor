function [datastructure,errorMsg] = SBimportCSVdata(filename)
% SBimportCSVdata
% Imports experimental data stored in an CSV (comma separated value) file
% into the datastructure used by the SBdata object.
% Please note that a special format of the CSV data is required. 
% This format is explained in the user's reference manual and example 
% files can be found in the SBTOOLBOX/examples folder.
% 
% filename: name of the .csv file containing the data
%
% datastructure: data structure used by SBdata object 
%                (empty if error occurred)
% errorMsg: string containing possible error messages. 

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
% INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
errorMsg = '';
datastructure = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize name and notes with default values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% use the filename as name for the data object
datastructure.name = strrep(filename,'.csv','');
datastructure.notes = 'no notes';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Open the file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(filename,'r');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for name
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The name is expected in the first row - first entry in this row needs to
% be 'Name:'
% get first row
dataRow = fgetl(fid);
% check if name identifier present
if ~isempty(strfind(lower(dataRow),'name:')),
    % name definition is present
    nameRow = strtrim(dataRow);
    name = strtrim(nameRow(6:end));
    if ~isempty(name),
        datastructure.name = name;
    end
    % name has been found, get new row from file to search for notes
    dataRow = fgetl(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for notes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The notes are expected in consecutvie rows, all starting with 'Notes:'
% check if notes identifier present
notes = '';
while ~isempty(strfind(lower(dataRow),'notes:')),
    % notes definition is present in current data row
    notesRow = strtrim(dataRow);
    notes = sprintf('%s%s\n',notes,strtrim(notesRow(7:end)));
    dataRow = fgetl(fid);
end
if ~isempty(notes),
    datastructure.notes = notes;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize time vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datastructure.time = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for stimulus flag
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The stimuli flag definitions are expected after the notes definitions in
% one row separated by commas and prefixed by 'StimulusFlag:'
% check if stimulus flag identifier present
stimulusFlagData = [];
if ~isempty(strfind(lower(dataRow),'stimulusflag:')),
    % stimulus flag definition is present
    stimulusFlagRow = strtrim(dataRow);
    % take away leading and trailing white space - no extra ',' needed
    % since a numeric value for each components is expected if row is present
    stimulusFlag = strtrim(stimulusFlagRow(14:end));
    if ~isempty(stimulusFlag),
        stimulusFlagData = strread(stimulusFlag,'%d','delimiter',',');
    end
    % check the contents of stimulusFlag (set all that is not "1" to "0")
    stimulusFlagData(find(stimulusFlagData~=1)) = 0;
    % stimulusFlag has been found, get new row from file to search for
    % noise offset
    dataRow = fgetl(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for noise offset
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The noise offset is expected after the stimulus Flag definitions in one
% row separated by commas and prefixed by 'NoiseOffset:'
% check if noise offset identifier present
noiseOffsetData = [];
if ~isempty(strfind(lower(dataRow),'noiseoffset:')),
    % noise offset definition is present
    noiseOffsetRow = strtrim(dataRow);
    % take away leading and trailing white space - no extra ',' needed
    % since a numeric value for each components is expected if row is present
    noiseOffset = strtrim(noiseOffsetRow(13:end));
    if ~isempty(noiseOffset),
        noiseOffsetData = strread(noiseOffset,'%f','delimiter',',');
    end
    % noiseOffset has been found, get new row from file to search for
    % noise variance
    dataRow = fgetl(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check for noise variance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The noise variance is expected after the noise offset definitions in one
% row separated by commas and prefixed by 'NoiseVariance:'
% check if noise variance identifier present
noiseVarianceData = [];
if ~isempty(strfind(lower(dataRow),'noisevariance:')),
    % noise variance definition is present
    noiseVarianceRow = strtrim(dataRow);
    % take away leading and trailing white space - no extra ',' needed
    % since a numeric value for each components is expected if row is present
    noiseVariance = strtrim(noiseVarianceRow(15:end));
    if ~isempty(noiseVariance),
        noiseVarianceData = strread(noiseVariance,'%f','delimiter',',');
    end
    % noiseVariance has been found, get new row from file to search for
    % component names
    dataRow = fgetl(fid);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Read the component names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The component names are expected in one single row, separated by commas
componentNames = strread(dataRow,'%s','delimiter',',');
if length(componentNames) == 0,
    errorMsg = sprintf('%s\nSyntax in component row seems to be wrong.\nIt should contain the names of the measured components, separated by commas\nand no identifier is required.\n',errorMsg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check length of noiseOffsetData against component names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(noiseOffsetData),
    if length(noiseOffsetData) ~= length(componentNames),
        errorMsg = sprintf('%s\nNoise offset data have been defined but their number does not match the number of components.\n',errorMsg);
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check length of noiseVarianceData against component names
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(noiseVarianceData),
    if length(noiseVarianceData) ~= length(componentNames),
        errorMsg = sprintf('%s\nNoise variance data have been defined but their number does not match the number of components.\n',errorMsg);
        return
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update component names and set other fields to default values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:length(componentNames),
    datastructure.measurements(k).name = componentNames{k};
    if ~isempty(stimulusFlagData),
        datastructure.measurements(k).stimulusflag = stimulusFlagData(k);
    else 
        datastructure.measurements(k).stimulusflag = 0;
    end        
    datastructure.measurements(k).data = [];
    if ~isempty(noiseOffsetData),
        datastructure.measurements(k).noiseoffset = noiseOffsetData(k);
    else
        datastructure.measurements(k).noiseoffset = 0;
    end
    if ~isempty(noiseVarianceData),
        datastructure.measurements(k).noisevariance = noiseVarianceData(k);
    else
        datastructure.measurements(k).noisevariance = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Import the time and experimental data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
experimentData = [];
while ~feof(fid),
    % read next data row
    dataRow = fgetl(fid);
    % add an extra comma at the end of the string
    dataRow = strcat(dataRow,',');
    % exchange all ',,' against ',NaN,'
    dataRow = strrep(dataRow,',,','NaN');
    % parse data
    dataRow = strread(dataRow,'%f','delimiter',',');
    % copy data into dataMatrix
    experimentData = [experimentData; [dataRow(1:end)]'];
end
% close the file
fclose(fid);
% extract the time and measurement data
timeData = experimentData(:,1);
% check the time vector
if ~isempty(find(isnan(timeData))),
    errorMsg = sprintf('%s\nNot all time instants defined in the time vector.\n',errorMsg);
end
test = timeData(2:end)-timeData(1:end-1);
if ~isempty(find(test <= 0)),
    errorMsg = sprintf('%s\nTime vector is not monotonically increasing.\n',errorMsg);
end
% add time vector to structure
datastructure.time = timeData;
% get measurement data
valueData = experimentData(:,2:length(componentNames)+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process the time and measurements
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for k = 1:length(componentNames),
    % extract the data from the cell-array-matrix
    datastructure.measurements(k).data = valueData(:,k);
end
return
