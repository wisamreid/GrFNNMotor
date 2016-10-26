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
% CONVERT EVENT DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% at this point the output data of this function need to be defined.
% the trigger expressions need to be converted to continuous functions.
% a trigger going from false to true is modelled here as a
% continuous function going from negative to positive.
% cycle through all events
valueString = '';
for k = 1:length(sbm.events),
    trigger = sbm.events(k).trigger;
    % check which operator used (le,lt,ge,gt)
    foundL = 0;
    foundG = 0;
    if ~isempty(strfind(trigger,'le')) || ~isempty(strfind(trigger,'lt')),
        foundL = 1;
    end
    if ~isempty(strfind(trigger,'ge')) || ~isempty(strfind(trigger,'gt')),
        foundG = 1;
    end
    if (foundL == 0 && foundG == 0) || (foundL == 1 && foundG == 1),
        error('The trigger expression for event ''%s'' is not correct.',sbm.events(k).name);
    end
    % get the expressions
    arguments = extractPSB(trigger);
    arguments = explodePCSB(arguments);
    if length(arguments) ~= 2,
        error('The trigger expression for event ''%s'' is not correct.',sbm.events(k).name);
    end
    % write the trigger string
    if foundL == 1,
        triggerString = sprintf('%s-%s',arguments{2},arguments{1});
    elseif foundG == 1,
        triggerString = sprintf('%s-%s',arguments{1},arguments{2});
    end
    % write the value string
    valueString = sprintf('%s%s\n',valueString,triggerString);
end
valueString = valueString(1:end-1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPEN FILE FOR WRITING AND WRITE HEADER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(filename,'w');
fprintf(fid,'function [value,isterminal,direction] = %s(time,statevector)\n',functionName);
fprintf(fid,'%% This function is passed to the MATLAB integrator in order to detect\n%% events that are present in the model.\n\n');

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
% WRITE THE EVENT DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'\n%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'%% EVENT DATA\n');
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'value = [%s]'';\n',valueString);
fprintf(fid,'isterminal = ones(1,%d);\n',length(sbm.events));
fprintf(fid,'direction = ones(1,%d);\n',length(sbm.events));
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