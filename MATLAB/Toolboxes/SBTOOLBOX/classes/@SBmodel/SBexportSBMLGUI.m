function SBexportSBMLGUI(model)
% SBexportSBMLGUI
% This function is more or less the counter part to "SBexportSBML" with the
% difference that ths function give you the opporunity to make last changes
% wihtin a Java GUI
%
% USAGE:
% ======
% SBexportSBMLGUI(model)
%
% model: SBmodel to export to SBML
% 
% changes made within the GUI can be revoked by pressing the "Cancel"
% button or they can be assigned to the variable "changedModel" available
% in the MATLAB environment after clicking the "Exit" button

% Information:
% ============
% Systems Biology Toolbox for MATLAB
% Copyright (C) 2005-2007 Fraunhofer-Chalmers Research Centre
% Author of the toolbox: Henning Schmidt, henning@hschmidt.de
% Programmed by Gunnar Drews, Student at University of
% Rostock Department of Computerscience, gunnar.drews@uni-rostock.de
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK IF SBML Toolbox is present
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~exist('OutputSBML'),
    error('The SBML Toolbox is not installed on your system. You need it for the SBML export.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK IF MODEL CONTAINS MATLAB FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(model.functionsMATLAB),
    error('Model contains MATLAB functions. Such functions are not supported in SBML.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test wether all name fields are filled and are unique
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
namesOK = checkNamesInSBmodel(model, 1);
if ~namesOK,
    errorMessage{1} = 'SBML export aborted because of errors concerning the names used in the given SBmodel. For further information use "checkNamesInSBmodel" function!';
    messageOutput(errorMessage, 1);
    return    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PREPROCESS and CHECK SBmodel (fill type, compartment, unittype fields if possible)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check wether model is a Java SBmodel or a MATLAB SBmodel
model = SBmakeSBMLpresets(model);
testSBMLsettings(model, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert Model to Java Object and start GUI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
javamodel = modelSB2Java(model);
SBexportSBMLGUI.SB2SBMLGUI(javamodel);

end