function addJars2DynamicCPsecure
% addStaticCP
% 
% This script utilizes setting the dynamic class path for the
% SBexportSBMLGUI and the therefor needed javamodel.
%
%
% USAGE
% ======
% 
% addJars2DynamicCP
%

% Information:
% ============
% Systems Biology Toolbox for MATLAB
% Copyright (C) 2005-2007 Fraunhofer-Chalmers Research Centre
% Author of the toolbox: Henning Schmidt, henning@hschmidt.de
% programmed by Gunnar Drews student at Depertment of Computerscience
% University of Rostock
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
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
% USA.

% define all known parts of paths
if isunix,
    locationOfJarsInSBTOOLBOX = '/gui/SBexportSBML/';
else
    locationOfJarsInSBTOOLBOX = '\gui\SBexportSBML\';
end
path2SBTOOLBOX = fileparts(which('installSB'));

javamodel = char([double(path2SBTOOLBOX), double(locationOfJarsInSBTOOLBOX), double('javamodel.jar')]);
SBexportSBMLGUI = char([double(path2SBTOOLBOX), double(locationOfJarsInSBTOOLBOX), double('SBexportSBMLGUI.jar')]);

isSetStaticModelJar = false;
isSetStaticGuiJar = false;
isSetDynamicModelJar = false;
isSetDynamicGuiJar = false;

% check wether one or both of our jars are already contained in static
% class path
tlineCounter = 1;
cpFileString = char([double(matlabroot), double('\toolbox\local\classpath.txt')]);
cptxtID = fopen(cpFileString, 'r');
while 1
    % read "classpath.txt" line per line
    tline = fgetl(cptxtID);
    % test for end of file
    if ~ischar(tline), break, end
    if length(strfind(tline, javamodel)),
        isSetStaticModelJar = true;
    end
    if length(strfind(tline, SBexportSBMLGUI)),
        isSetStaticGuiJar = true;
    end
end
fclose(cptxtID);

% if jars not have been declared static look in dynamic path
if (~isSetStaticModelJar || ~isSetStaticGuiJar),
    dynPaths = javaclasspath('-dynamic');
    for index = 1 : length(dynPaths),
        if length(strfind(dynPaths{index}, javamodel)),
            isSetDynamicModelJar = true;
        end
        if length(strfind(dynPaths{index}, SBexportSBMLGUI)),
            isSetDynamicGuiJar = true;
        end
    end
end

% if one or both of our jars haven't been added to path
% add now to dynamic java class path
if (~isSetStaticModelJar && ~isSetDynamicModelJar),
    javaaddpath(javamodel);
end
if (~isSetStaticGuiJar && ~isSetDynamicGuiJar),
    javaaddpath(SBexportSBMLGUI);
end

end % of addJars2DynamicCP