function [SBstructure,errorMsg] = SBimportSBML(SBMLmodelFile)
% SBimportSBML 
% imports a SBML model using the TranslateSBML function from libSBML
% Supported SBML levels: 1 and 2
%
% SBMLmodelFile: SBMLmodelfilename.xml (string)
%
% SBstructure: empty if error occurred

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
% INITIALIZATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
errorMsg = '';
SBstructure = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SBML -> MATLAB
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use the TranslateSBML function from libSBML
SBMLmodel = TranslateSBML(SBMLmodelFile);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONVERT INTO OWN SB STRUCTURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp(SBMLmodel.typecode,'SBML_MODEL'),
    errorMsg = 'Model is not a SBML model';
    return;
end
if SBMLmodel.SBML_level == 1,
    % Convert level 1
    [SBstructure,errorMsg] = SBconvertSBML1(SBMLmodel);
elseif SBMLmodel.SBML_level == 2 && SBMLmodel.SBML_version == 1,
    % Convert level 2
    [SBstructure,errorMsg] = SBconvertSBML2(SBMLmodel);
else
    % Not supported SBML level
    errorMsg = sprintf('Level %d Version %d of SBML not yet implemented',SBMLmodel.SBML_level,SBMLmodel.SBML_version);
end
