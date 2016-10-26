function [] = SBexportXLSdata(varargin)
% SBexportXLSdata
% Exports an SBdata object to a XLS (excel) file.
% The format of the written XLS file is explained in the user's reference
% manual and example files can be found in the SBTOOLBOX/examples folder.
% 
% USAGE:
% ======
% [] = SBexportXLSdata(data)
% [] = SBexportXLSdata(data,filename)
%
% data: SBdata object containing the data
% filename: desired filename for XLS file. The extension '.xls' is not
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
filename = strcat(filename,'.xls');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK IF OBJECT CONTAINS MEASUREMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numberComponents = length(data.measurements);
if numberComponents == 0,
    error('The object does not contain any measurements');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Construct Cell-Array Matrix of correct size
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The number of columns is numberComponents+1
% The number of rows is 6 (rows header) + length of time vector
RAW = cell(6+length(data.time),numberComponents+1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write identifiers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RAW{1,1} = 'Name';
RAW{2,1} = 'Notes';
RAW{3,1} = 'Components';
RAW{4,1} = 'StimulusFlag';
RAW{5,1} = 'NoiseOffset';
RAW{6,1} = 'NoiseVariance';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fill Matrix with information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First the header data
RAW{1,2} = data.name;
RAW{2,2} = data.notes;
RAW(7:7+length(data.time)-1,1) = num2cell(data.time);
for k = 1:numberComponents,
    RAW{3,k+1} = data.measurements(k).name;
    RAW{4,k+1} = data.measurements(k).stimulusflag;
    RAW{5,k+1} = data.measurements(k).noiseoffset;
    RAW{6,k+1} = data.measurements(k).noisevariance;
    RAW(7:7+length(data.measurements(k).data)-1,k+1) = num2cell(data.measurements(k).data);
end
[success,message] = xlswrite(filename,RAW);
return
