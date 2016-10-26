function writeSBML(model)
% writeSBML
%
% This function serves as shortcut for writing a model into a SBML file
% It's mainly used by the SBexportSBMLGUI which call this script from Java.
%
% USAGE:
% ======
% writeSBML(model)
%
% model: SBmodel or SBmodelJava that will be converted into a MATLAB-SBML
%        structure and is written into a SBML file afterwards

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

% check wether model is a Java SBmodel or a MATLAB SBmodel
if isjava(model),
    SBmodel = java2modelSB(model);
else
    SBmodel = model;
end

% use export script to convert SBmodel to MATLAB SBML structure and write
% it into SBML file
SBexportSBML(SBmodel, 1);

end