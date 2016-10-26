function addJars2DynamicCP
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

% suppress warnings for already sttic defined class paths
warning off MATLAB:javaclasspath:jarAlreadySpecified;

% add Jars to dynamic
javaaddpath(javamodel);
javaaddpath(SBexportSBMLGUI);

% turn warnings on again
warning on MATLAB:javaclasspath:jarAlreadySpecified;

end % of addJars2DynamicCP