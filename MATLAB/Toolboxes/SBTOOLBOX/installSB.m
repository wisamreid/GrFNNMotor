% installSB
% Installation script for the SBtoolbox. Edit the data below to match your
% system and run it.

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

% clean up workspace and desktop
clear classes; clc; close all;
VERSION = 'Version 1.8';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EDIT THE FOLLOWING VARIABLES TO MATCH YOUR SYSTEM
% Alternatively you can add these two paths manually!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PATH_SBMLTOOLBOX = 'C:\Program Files\SBML\SBMLToolbox-2.0.2';
PATH_XPPAUT = 'C:\xppall';

% apply path settings 
PATH_SBTOOLBOX = fileparts(which('installSB.m'));
addpath(genpath(PATH_SBTOOLBOX))
addpath(genpath(PATH_SBMLTOOLBOX))
addpath(genpath(PATH_XPPAUT))
% Add path to temporary directory
addpath(tempdir)
% save the path changes
savepath

% add SBML export jars to javaclasspath (static if windows, dynamic if unix)
if isunix,
    addJars2DynamicCP
else
    addJars2StaticCP
end

% output license information.
disp(sprintf('Systems Biology Toolbox for MATLAB - %s\nCopyright (C) 2005-2006 Henning Schmidt, FCC, henning@fcc.chalmers.se\nThis is free software, distributed under the GNU General Public License.\nYou are welcome to redistribute it under certain conditions; type\n''open copying.html'' for details. The toolbox comes with ABSOLUTELY NO\nWARRANTY; for details type ''open warranty.html''. The complete GNU GPL\ncan be found at http://www.gnu.org/licenses/gpl.html.\nThe newest version can always be downloaded from http://www.sbtoolbox.org\n',VERSION));