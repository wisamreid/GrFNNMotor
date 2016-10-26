function [datastructure,errorMsg] = SBimportXLSdata(filename)
% SBimportXLSdata
% Imports experimental data stored in an XLS Excel file
% into the datastructure used by the SBdata object.
% Please note that a special format of the excel data is required. 
% This format is explained in the user's reference manual and example 
% files can be found in the SBTOOLBOX/examples folder.
% 
% filename: name of the .xls file containing the data
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
% Import the Excel data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[NUMERIC,TXT,RAW] = xlsread(filename);
% Skip the NUMERIC and TEXT variables and only consider the RAW data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process identifier names and check errors
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% - each identifier can only appear once
% - Components and Measurements identifier are required
% - Name, Notes, Formula, NoiseOffset, and NoiseVariance are optional
foundName = 0;
foundNotes = 0;
foundComponents = 0;
foundStimulusFlag = 0;
foundNoiseOffset = 0;
foundNoiseVariance = 0;
foundMeasurements = 0;
k = 1; 
while ~isnumeric(RAW{k,1}),
    % search for name identifier in first column
    if strcmp('name',lower(removeWhiteSpace(RAW{k,1})))
        if foundName == 1,
            errorMsg = sprintf('%s\nThe ''Name'' identifier has been found more than once.\n',errorMsg);
        end
        % get row number for name
        rowName = k;
        foundName = 1;
    end
    % search for notes identifier in first column
    if strcmp('notes',lower(removeWhiteSpace(RAW{k,1})))
        if foundNotes == 1,
            errorMsg = sprintf('%s\nThe ''Notes'' identifier has been found more than once.\n',errorMsg);
        end
        % get row number for notes
        rowNotes = k;
        foundNotes = 1;
    end
    % search for Components identifier in first column
    if strcmp('components',lower(removeWhiteSpace(RAW{k,1})))
        if foundComponents == 1,
            errorMsg = sprintf('%s\nThe ''Components'' identifier has been found more than once.\n',errorMsg);
        end
        % get row number for Components
        rowComponents = k;
        foundComponents = 1;
    end
    % search for Components identifier in first column
    if strcmp('stimulusflag',lower(removeWhiteSpace(RAW{k,1})))
        if foundStimulusFlag == 1,
            errorMsg = sprintf('%s\nThe ''StimulusFlag'' identifier has been found more than once.\n',errorMsg);
        end
        % get row number for Components
        rowStimulusFlag = k;
        foundStimulusFlag = 1;
    end
    % search for NoiseOffset identifier in first column
    if strcmp('noiseoffset',lower(removeWhiteSpace(RAW{k,1})))
        if foundNoiseOffset == 1,
            errorMsg = sprintf('%s\nThe ''NoiseOffset'' identifier has been found more than once.\n',errorMsg);
        end
        % get row number for NoiseOffset
        rowNoiseOffset = k;
        foundNoiseOffset = 1;
    end
    % search for NoiseVariance identifier in first column
    if strcmp('noisevariance',lower(removeWhiteSpace(RAW{k,1})))
        if foundNoiseVariance == 1,
            errorMsg = sprintf('%s\nThe ''NoiseVariance'' identifier has been found more than once.\n',errorMsg);
        end
        % get row number for NoiseVariance
        rowNoiseVariance = k;
        foundNoiseVariance = 1;
    end
    k = k+1;
end
if foundComponents == 0,
    errorMsg = sprintf('%s\nThe ''Components'' identifier has not been found.\n',errorMsg);
end
if foundNoiseOffset ~= foundNoiseVariance,
    errorMsg = sprintf('%s\nIf ''NoiseOffset'' is defined, also ''NoiseVariance'' has to be defined, and vice-versa.\n',errorMsg);
end
startRowData = k;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% If error occurred so far return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(errorMsg),
    return;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process the name (Name is OPTIONAL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the name is supposed to be in the second column
if foundName,
    datastructure.name = RAW{rowName,2};
else
    datastructure.name = 'untitled';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process the notes (Notes is OPTIONAL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% the notes data is supposed to be in the second column
if foundNotes,
    datastructure.notes = RAW{rowNotes,2};
else
    datastructure.notes = 'no notes';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize time field in structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
datastructure.time = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process the component names (Components are REQUIRED)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% insert component names in datastructure (skip NANs)
% the number of component names then determines the number 
% of measurements that are assumed to be present in the file
numberComponents = 0;
while (numberComponents+2 <= size(RAW,2)),
    datastructure.measurements(numberComponents+1).name = RAW{rowComponents,numberComponents+2};
    numberComponents = numberComponents + 1;
end
% check that data is present
if numberComponents == 0,
    errorMsg = sprintf('%s\nNo components present in the data.\n',errorMsg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process the stimulusFlag information (optional)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% stimulus flag information is optional - if given then assume that all 
% elements are measurements (value = 0)
for k = 1:numberComponents,
    if foundStimulusFlag && RAW{rowStimulusFlag,k+1} == 1,
        datastructure.measurements(k).stimulusflag = 1;
    else
        datastructure.measurements(k).stimulusflag = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process the noise information (NoiseOffset and NoiseVariance are OPTIONAL)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% noise information is optional - if something is not given then use a
% default value of 0. The noise data needs to be specified in the same
% column as the component names.
for k = 1:numberComponents,
    if foundNoiseOffset && sum(~isnan(RAW{rowNoiseOffset,k+1})),
        datastructure.measurements(k).noiseoffset = RAW{rowNoiseOffset,k+1};
    else
        datastructure.measurements(k).noiseoffset = 0;
    end
    if foundNoiseVariance && sum(~isnan(RAW{rowNoiseVariance,k+1})),
        datastructure.measurements(k).noisevariance = RAW{rowNoiseVariance,k+1};
    else
        datastructure.measurements(k).noisevariance = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process the measurement information (Measurement data is REQUIRED)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
timeData = RAW(startRowData:end,1);
timeData = [timeData{:}]';
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
% get the measurement data
experimentData = RAW(startRowData:end,2:numberComponents+1);
% cycle through the measurementData and extract data to put into the
% datastructure
for k = 1:numberComponents,
    % extract the data from the cell-array-matrix
    valueData = experimentData(:,k);
    % convert data from cell-array to vector
    valueVector = [valueData{:}]';
    % copy data into structure
    datastructure.measurements(k).data = valueVector;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DELETE WHITESPACES IN STRINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Useful for taking away whitespaces in kineticLaw formulas, as
% seen in some example models
function [outputString] = removeWhiteSpace(inputString)
outputString = strrep(inputString,' ','');
return

