function [] = display(sbm)
% display: Displays information about SBmodel. This function is 
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
% COLLECT INFORMATION ABOUT THE MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numberStates = length(sbm.states);
numberVariables = length(sbm.variables);
numberParameters = length(sbm.parameters);
numberReactions = length(sbm.reactions);
numberFunctions = length(sbm.functions);
numberEvents = length(sbm.events);
functionsMATLABpresent = ~isempty(sbm.functionsMATLAB);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY INFORMATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
text = sprintf('\tSBmodel\n\t=======\n');
text = sprintf('%s\tName: %s\n',text,sbm.name);
text = sprintf('%s\tNumber States:\t\t%d\n',text,numberStates);
text = sprintf('%s\tNumber Variables:\t%d\n',text,numberVariables);
text = sprintf('%s\tNumber Parameters:\t%d\n',text,numberParameters);
text = sprintf('%s\tNumber Reactions:\t%d\n',text,numberReactions);
text = sprintf('%s\tNumber Functions:\t%d\n',text,numberFunctions);
if numberEvents > 0,
    text = sprintf('%s\tNumber Events:\t%d\n',text,numberEvents);
end
if functionsMATLABpresent,
    text = sprintf('%s\tMATLAB functions present',text);
end
disp(text);