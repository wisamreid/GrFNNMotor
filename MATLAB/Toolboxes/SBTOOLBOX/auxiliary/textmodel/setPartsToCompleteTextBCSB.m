function [completeTextBC] = setPartsToCompleteTextBCSB(modelTextStructureBC)
% setPartsToCompleteTextBCSB: Sets the different parts of a biochemical
% oriented text description of an SBmodel together to the complete text

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
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
% USA.

completeTextBC = sprintf('********** MODEL NAME\n%s\n',modelTextStructureBC.name);
completeTextBC = sprintf('%s\n********** MODEL NOTES\n%s\n',completeTextBC,modelTextStructureBC.notes);
completeTextBC = sprintf('%s\n********** MODEL STATE INFORMATION\n%s\n',completeTextBC,modelTextStructureBC.states);
completeTextBC = sprintf('%s\n********** MODEL PARAMETERS\n%s\n',completeTextBC,modelTextStructureBC.parameters);
completeTextBC = sprintf('%s\n********** MODEL VARIABLES\n%s\n',completeTextBC,modelTextStructureBC.variables);
completeTextBC = sprintf('%s\n********** MODEL REACTIONS\n%s\n',completeTextBC,modelTextStructureBC.reactions);
completeTextBC = sprintf('%s\n********** MODEL FUNCTIONS\n%s\n',completeTextBC,modelTextStructureBC.functions);
completeTextBC = sprintf('%s\n********** MODEL EVENTS\n%s\n',completeTextBC,modelTextStructureBC.events);
completeTextBC = sprintf('%s\n********** MODEL MATLAB FUNCTIONS\n%s\n',completeTextBC,modelTextStructureBC.functionsMATLAB);
return
