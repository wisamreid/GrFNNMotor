function [] = display(data)
% display: Displays information about SBdata. This function is 
% called by MATLAB whenever an object is the result of a statement that
% is not terminated by a semicolon. 

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
% COLLECT INFORMATION ABOUT THE DATA OBJECT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
name = data.name;
notes = data.notes;
time = data.time;
numberComponentsMeasured = 0;
numberComponentsStimuli = 0;
noiseInformation = 0;
allMeasurementsDone = 1;
for k = 1:length(data.measurements),
    if data.measurements(k).stimulusflag == 0,
        numberComponentsMeasured = numberComponentsMeasured + 1;
    else
        numberComponentsStimuli = numberComponentsStimuli + 1;
    end
    % noise information
    if data.measurements(k).noiseoffset ~= 0 || data.measurements(k).noisevariance ~= 0,
        noiseInformation = noiseInformation + 1;
    end
    % check if NaN present in data
    if ~isempty(find(isnan(data.measurements(k).data))),
        allMeasurementsDone = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY INFORMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
text = sprintf('\tSBdata\n\t======\n');
text = sprintf('%s\tName: %s\n',text,name);
text = sprintf('%s\tMeasured components:\t%d\n',text,numberComponentsMeasured);
text = sprintf('%s\tStimuli components:\t\t%d\n',text,numberComponentsStimuli);
text = sprintf('%s\tNumber time points: \t%d\n',text,length(time));
if noiseInformation == 0,
    text = sprintf('%s\tNo noise information present\n',text);
else
    text = sprintf('%s\tNoise information present for %d components\n',text,noiseInformation);
end
if allMeasurementsDone == 1,
    text = sprintf('%s\tMeasurements for all time points present\n',text);
else
    text = sprintf('%s\tMeasurements not present for all time points\n',text);
end
disp(text);