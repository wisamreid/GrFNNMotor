 function [] = SBcreateEventFunction(sbm,filename)
% SBcreateEventFunction: writes the event handling function that is
% necessary to pass to the intergrator in order to be able to deal with
% discrete state events when doing simulations. This function is called is
% called by the function SBcreateODEfile in the case that the event flag is
% set and at least one event is present in the model.
%
% USAGE:
% ======
% [] = SBcreateEventFunction(sbm,filename)
%
% sbm: SBmodel  (ODE file model description is not possible to use)
% filename: the filename to which wo write the model

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
% CHECK IF SBmodel (not really needed but safer)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~strcmp('SBmodel',class(sbm)),
    error('Function only defined for SBmodels.');
end

[PATHSTR,functionName,EXT,VERSN] = fileparts(filename);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPEN FILE FOR WRITING AND WRITE HEADER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(filename,'w');
fprintf(fid,'function [newstates] = %s(eventIndex,time,statevector)\n',functionName);
fprintf(fid,'%% This function is used to determine the resulting new states when an event has fired.\n\n');

fprintf(fid,'newParameters = 0;\n');
fprintf(fid,'global modelParameterValuesNew\n');
fprintf(fid,'if ~isempty(modelParameterValuesNew),\n');
fprintf(fid,'    parameterValuesNew = modelParameterValuesNew;\n');
fprintf(fid,'    newParameters = 1;\n');
fprintf(fid,'end\n\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITE THE MODEL STUFF TO INITIALIZE all eventually needed data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'%% MODEL DATA\n');
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
% PROCESS STATES
stateNames = SBstates(sbm);
for k = 1:length(stateNames),
   fprintf(fid,'%s = statevector(%d);\n',stateNames{k},k);
end
% PROCESS PARAMETERS
[parameterNames,parameterValues] = SBparameters(sbm);
fprintf(fid,'if newParameters == 0,\n');
for k = 1:length(parameterNames),
    fprintf(fid,'\t%s = %g;\n',parameterNames{k},parameterValues(k));
end
fprintf(fid,'else\n'); 
for k = 1:length(parameterNames),
    fprintf(fid,'\t%s = parameterValuesNew(%d);\n',parameterNames{k},k);
end
fprintf(fid,'end\n');
% PROCESS VARIABLES
[variableNames,variableFormulas] = SBvariables(sbm);
for k = 1:length(variableNames),
   fprintf(fid,'%s = %s;\n',variableNames{k},variableFormulas{k});
end 
% PROCESS REACTIONS
[reactionNames,reactionFormulas] = SBreactions(sbm);
for k = 1:length(reactionNames),
   fprintf(fid,'%s = %s;\n',reactionNames{k},reactionFormulas{k});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXECUTE THE EVENT ASSIGNMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'%% EVENT ASSIGNMENTS\n');
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
% assign new states with old states
for k = 1:length(stateNames),
   fprintf(fid,'%s_new = %s;\n',stateNames{k},stateNames{k});
end

% then do the event assignments
for k = 1:length(sbm.events),
    fprintf(fid,'if sum(ismember(eventIndex,%d)) ~= 0,\n',k);
    for k2 = 1:length(sbm.events(k).assignment),
        fprintf(fid,'%s_new = %s;\n',sbm.events(k).assignment(k2).variable,sbm.events(k).assignment(k2).formula);
    end
    fprintf(fid,'end\n');
end

% then return the changed new states
for k = 1:length(stateNames),
   fprintf(fid,'newstates(%d) = %s_new;\n',k,stateNames{k});
end
fprintf(fid,'\n');

% return
fprintf(fid,'return\n');
fprintf(fid,'\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Continue writing model stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITE MODEL FUNCTIONS
[functionNames,functionFormulas,functionArguments] = SBfunctions(sbm);
for k = 1:length(functionNames),
    fprintf(fid,'function [result] = %s(%s)\n',functionNames{k},functionArguments{k});
    fprintf(fid,'result = %s;\n',functionFormulas{k});
    fprintf(fid,'return\n');
    fprintf(fid,'\n');
end
% WRITE THE MATLAB FUNCTIONS
functionsMATLAB = SBfunctionsMATLAB(sbm);
if ~isempty(functionsMATLAB),
    fprintf(fid,'%s',functionsMATLAB);
    fprintf(fid,'\n');
end
% CLOSE FILE
fclose(fid);
return