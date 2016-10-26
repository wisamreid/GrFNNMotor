function [] = SBexportCSVdata(varargin)
% SBexportCSVdata
% Exports an SBdata object to a CSV (comma separated value) file.
% The format of the written CSV file is explained in the user's reference
% manual and example files can be found in the SBTOOLBOX/examples folder.
% 
% USAGE:
% ======
% [] = SBexportCSVdata(data)
% [] = SBexportCSVdata(data,filename)
%
% data: SBdata object containing the data
% filename: desired filename for CSV file. The extension '.csv' is not
%   required.
%
% DEFAULT VALUES:
% ===============
% filename: constructed from the data objects name

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
    data = varargin{1};
    % if no filename provided then use the name of the SBdata object 
    % as filename 
    filename = data.name;
elseif nargin == 2,
    data = varargin{1};
    % extract filename from input arguments to skip eventual extension
    [PATHSTR,filename,EXT,VERSN] = fileparts(varargin{2});
else
    error('Incorrect number of input arguments.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS FILENAME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
filename = strrep(filename,' ','');
filename = strrep(filename,'-','');
filename = strrep(filename,'+','');
filename = strrep(filename,'*','');
filename = strrep(filename,'/','');
filename = strcat(filename,'.csv');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK IF OBJECT CONTAINS MEASUREMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(data.measurements) == 0,
    error('The object does not contain any measurements');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct Name String
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nameString = sprintf('Name: %s',data.name);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct Notes String
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
notesString = strtrim(sprintf('Notes: %s',data.notes));
notesString = strrep(notesString,sprintf('\n'),sprintf('\nNotes: '));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct StimulusFlag, Component, NoiseOffset, and NoiseVariance String
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numberComponents = length(data.measurements);
componentNames = {};
stimulusFlagString = 'StimulusFlag: ';
componentString = '';
noiseOffsetString = 'NoiseOffset: ';
noiseVarianceString = 'NoiseVariance: ';
maxDataLength = 0;
for k = 1:numberComponents,
    componentNames{k} = data.measurements(k).name;
    if isempty(componentNames{k}),
        error('The %d-th component name is undefined.');
    end
    stimulusFlagString = sprintf('%s%d,',stimulusFlagString,data.measurements(k).stimulusflag);
    componentString = sprintf('%s%s,',componentString,componentNames{k});
    noiseOffsetString = sprintf('%s%f,',noiseOffsetString,data.measurements(k).noiseoffset);
    noiseVarianceString = sprintf('%s%f,',noiseVarianceString,data.measurements(k).noisevariance);
end
% take away the last comma in some strings
stimulusFlagString = stimulusFlagString(1:end-1);
componentString = componentString(1:end-1);
noiseOffsetString = noiseOffsetString(1:end-1);
noiseVarianceString = noiseVarianceString(1:end-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create the file and write the information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(filename,'w');
fprintf(fid,'%s\n',nameString);
fprintf(fid,'%s\n',notesString);
fprintf(fid,'%s\n',stimulusFlagString);
fprintf(fid,'%s\n',noiseOffsetString);
fprintf(fid,'%s\n',noiseVarianceString);
fprintf(fid,'%s\n',componentString);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct the data matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
experimentData = data.time(:);
for k = 1:numberComponents,
    experimentData(:,k+1) = data.measurements(k).data(:);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write data to string and close file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataString = '';
for k = 1:size(experimentData,1),
    dataString = sprintf('%f,',experimentData(k,:));
    fprintf(fid,'%s\n',dataString(1:end-1));
end
fclose(fid);
return
