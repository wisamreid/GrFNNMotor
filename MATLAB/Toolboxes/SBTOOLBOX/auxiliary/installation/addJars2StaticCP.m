function addJars2StaticCP
% addStaticCP
% 
% This script utilizes setting the static class path for the
% SBexportSBMLGUI and the therefor needed javamodel.
%
%
% USAGE
% ======
% 
% addJars2StaticCP
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

javamodel = 'javamodel.jar';
SBexportSBMLGUI = 'SBexportSBMLGUI.jar';
% define all known parts of paths
if isunix,
    locationOfJarsInSBTOOLBOX = '/gui/SBexportSBML/';
else
    locationOfJarsInSBTOOLBOX = '\gui\SBexportSBML\';
end
path2SBTOOLBOX = fileparts(which('installSB'));

% define end of line character(s) for rebuilding "classpath.txt"
% depending on used machine
if ispc, eoL = sprintf('\r\n'); end
if isunix, eoL = sprintf('\n'); end

jarfound = false;

tlineCounter = 1;
cpFileString = char([double(matlabroot), double('\toolbox\local\classpath.txt')]);
cptxtID = fopen(cpFileString, 'r');
% jars may have already been added to static class path
% take a look
while 1
    % in the case jars are already contained with a
    % false path we have to rebuild "classpath.txt"
    % therefor collect all lines except the ones
    % containing our jars
    tline = fgets(cptxtID);
    % test for end of file
    if ~ischar(tline), break, end
    if (length(strfind(tline, javamodel)) || length(strfind(tline, SBexportSBMLGUI))),
        jarfound = true;
    elseif strcmp(eoL, tline),
        % skip empty lines
    else
        tlineArray{tlineCounter} = tline;
        tlineCounter = tlineCounter + 1;
    end
end
fclose(cptxtID);

if jarfound,
    % entries for jar files have been found
    % rewriting "classpaath.txt" seems to be easier
    % than testing wether the entries may be correct
    % and manipulating them bytewise
    
    % open "classpath.txt" for writing and delete former content
    cptxtID = fopen(cpFileString, 'w');
    for index = 1 : length(tlineArray),
        fwrite(cptxtID, tlineArray{index});
    end
    % now add jars
    fwrite(cptxtID, char([double(path2SBTOOLBOX), double(locationOfJarsInSBTOOLBOX), double(javamodel), double(eoL)]));
    fwrite(cptxtID, char([double(path2SBTOOLBOX), double(locationOfJarsInSBTOOLBOX), double(SBexportSBMLGUI), double(eoL)]));
    fclose(cptxtID);
else
    % the needed jars haven't been found in "classpath.txt",
    % it's sufficient to append them
    
    % open "classpath.txt" for appending
    cptxtID = fopen(cpFileString, 'a');
    fwrite(cptxtID, char(double(eoL)));
    fwrite(cptxtID, char([double(path2SBTOOLBOX), double(locationOfJarsInSBTOOLBOX), double(javamodel), double(eoL)]));
    fwrite(cptxtID, char([double(path2SBTOOLBOX), double(locationOfJarsInSBTOOLBOX), double(SBexportSBMLGUI), double(eoL)]));
    fclose(cptxtID);
end