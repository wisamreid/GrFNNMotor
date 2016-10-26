function applyChanges(sbmj)
% applyChanges
%
% This is only a helper script to make the changes made in the Java GUI
% available through a SBmodel
%
% USAGE:
% ======
% applyChanges(model)
%
% model: SBmodelJava that will be converted into a SBmodel and put into
%        variable for access after closing the SBexportSBML Java GUI

% Information:
% ============
% Systems Biology Toolbox for MATLAB
% Copyright (C) 2005-2007 Fraunhofer-Chalmers Research Centre
% Author of the toolbox: Henning Schmidt, henning@hschmidt.de
% Programmed by Gunnar Drews, gunnar.drews@uni-rostock.de
% Student at University of Rostock Department of Computerscience
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

global changedModel;

% check wether model is a Java SBmodel or a MATLAB SBmodel
if isjava(sbmj),
    changedModel = java2modelSB(sbmj);
else
    disp('Sorry, but this is only a helper script to make the changes made in the Java GUI available through a SBmodel!');
    return;
end

end