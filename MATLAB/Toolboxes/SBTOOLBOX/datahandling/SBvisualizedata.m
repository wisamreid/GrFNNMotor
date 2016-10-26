function [] = SBvisualizedata(data)
% SBvisualizedata
% Function allowing to visualize the content of an SBdata object. 
% The function just prepares the data. Display is then realized using
% the SBplot function.
% 
% USAGE:
% ======
% SBvisualizedata(data)
%
% data: SBdata object containing the data

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

% convert SBdata object to struct
data = SBstruct(data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK IF MODEL CONTAINS MEASUREMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(data.measurements) == 0,
    error('The model does not contain any measurements');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET TIME VECTOR
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
time = data.time;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS THE DATA INTO A MATRIX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% also construct dataNames and legendtext data
dataNames = {};
dataMatrix = NaN*ones(length(data.time),length(data.measurements));
timeMatrix = dataMatrix;
legendtext = {};
for k = 1:length(data.measurements),
    if data.measurements(k).stimulusflag == 1,
        dataNames{k} = sprintf('%s (stimulus)',data.measurements(k).name);
        legendtext{k} = sprintf('%s no/nv: %g/%g (stimulus)',data.measurements(k).name,data.measurements(k).noiseoffset,data.measurements(k).noisevariance);
    else
        dataNames{k} = data.measurements(k).name;
        legendtext{k} = sprintf('%s no/nv: %g/%g (measurement)',data.measurements(k).name,data.measurements(k).noiseoffset,data.measurements(k).noisevariance);
    end
    dataComponent = data.measurements(k).data;
    timeComponent = data.time;
    indicesMeasured = find(~isnan(dataComponent));
    dataComponent = dataComponent(indicesMeasured);
    timeComponent = timeComponent(indicesMeasured);
    dataMatrix(1:length(dataComponent),k) = dataComponent;
    timeMatrix(1:length(timeComponent),k) = timeComponent;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT DATA USING SBplot FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SBplot(timeMatrix,dataMatrix,dataNames,legendtext,'o-');
return