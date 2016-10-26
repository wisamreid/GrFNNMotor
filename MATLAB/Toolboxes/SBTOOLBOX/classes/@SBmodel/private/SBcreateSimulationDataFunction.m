function [] = SBcreateSimulationDataFunction(sbm,filename)
% SBcreateSimulationDataFunction: is called by SBcreateODEfile to create a
% function that is able to calculate the values of the variables and
% parameters for given state or state time series.
%
% USAGE:
% ======
% [] = SBcreateSimulationDataFunction(sbm, filename)
%
% sbm: SBmodel
% filename: the name of the file to be created

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
fprintf(fid,'function [output] = %s(varargin)\n',functionName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITE THE HEADER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'%% %s\n',sbm.name);
fprintf(fid,'%% Generated: %s\n',datestr(now));
fprintf(fid,'%% \n');
fprintf(fid,'%% [output] = %s(statevector) => time=0\n',functionName);
fprintf(fid,'%% [output] = %s(time,statevector)\n',functionName);
fprintf(fid,'%% \n');
fprintf(fid,'%% output: structure containing information about the values of the variables\n');
fprintf(fid,'%%         and reaction rates for given time and state information.\n');
fprintf(fid,'%% \n');
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLE VARARGINS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'%% HANDLE VARIABLE INPUT ARGUMENTS\n');
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'if nargin == 1,\n');

    fprintf(fid,'\ttime = [0];\n');
    fprintf(fid,'\tstatevector = varargin{1};\n');

fprintf(fid,'elseif nargin == 2,\n');

    fprintf(fid,'\ttime = varargin{1};\n');
    fprintf(fid,'\tstatevector = varargin{2};\n');

fprintf(fid,'else\n');
    fprintf(fid,'\terror(''Incorrect number of input arguments.'');\n');
fprintf(fid,'end\n');
fprintf(fid,'\n');

fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'%% DETERMINE DATA\n');
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'variableValuesTotal = [];\n');
fprintf(fid,'reactionValuesTotal = [];\n');
fprintf(fid,'if length(time) == 1,\n');
    fprintf(fid,'\t[variableValuesTotal, reactionValuesTotal] = getValues(time,statevector);\n');
fprintf(fid,'else\n');
    fprintf(fid,'\tfor k = 1:length(time),\n');
        fprintf(fid,'\t\t[variableValues, reactionValues] = getValues(time(k),statevector(k,:));\n');
        fprintf(fid,'\t\tvariableValuesTotal = [variableValuesTotal; variableValues];\n');
        fprintf(fid,'\t\treactionValuesTotal = [reactionValuesTotal; reactionValues];\n');
    fprintf(fid,'\tend\n');
fprintf(fid,'end\n');

variableNames = SBvariables(sbm);
reactionNames = SBreactions(sbm);
stateNames = SBstates(sbm);

fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'%% CONSTRUCT OUTPUT VARIABLE\n');
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'output.time = time;\n');
fprintf(fid,'output.states = ');
if length(stateNames) > 0,
    fprintf(fid,'{');
    for k = 1:length(stateNames),
        fprintf(fid,'''%s'',',stateNames{k});
    end
    fprintf(fid,'};\n');
else
    fprintf(fid,'{};\n');
end
fprintf(fid,'output.statevalues = statevector;\n');
fprintf(fid,'output.variables = ');
if length(variableNames) > 0,
    fprintf(fid,'{');
    for k = 1:length(variableNames),
        fprintf(fid,'''%s'',',variableNames{k});
    end
    fprintf(fid,'};\n');
else
    fprintf(fid,'{};\n');
end
fprintf(fid,'output.variablevalues = variableValuesTotal;\n');
fprintf(fid,'output.reactions = ');
if length(reactionNames) > 0,
    fprintf(fid,'{');
    for k = 1:length(reactionNames),
        fprintf(fid,'''%s'',',reactionNames{k});
    end
    fprintf(fid,'};\n');
else 
    fprintf(fid,'{};\n');
end
fprintf(fid,'output.reactionvalues = reactionValuesTotal;\n');

fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'%% RETURN\n');
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'return\n');
fprintf(fid,'\n');
fprintf(fid,'\n');

fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'%% getValues FUNCTION\n');
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'function [variableValues, reactionValues] = getValues(time,statevector)\n');
fprintf(fid,'variableValues = [];\n');
fprintf(fid,'reactionValues = [];\n');

% PROCESS STATES
fprintf(fid,'%% STATES\n');
stateNames = SBstates(sbm);
for k = 1:length(stateNames),
   fprintf(fid,'%s = statevector(%d);\n',stateNames{k},k);
end

% PROCESS PARAMETERS
fprintf(fid,'%% PARAMETERS\n');
[parameterNames,parameterValues] = SBparameters(sbm);
for k = 1:length(parameterNames),
    fprintf(fid,'%s = %g;\n',parameterNames{k},parameterValues(k));
end

% PROCESS VARIABLES
fprintf(fid,'%% VARIABLES\n');
[variableNames,variableFormulas] = SBvariables(sbm);
for k = 1:length(variableNames),
   fprintf(fid,'%s = %s;\n',variableNames{k},variableFormulas{k});
end

% PROCESS REACTIONS
fprintf(fid,'%% REACTIONS\n');
[reactionNames,reactionFormulas] = SBreactions(sbm);
for k = 1:length(reactionNames),
   fprintf(fid,'%s = %s;\n',reactionNames{k},reactionFormulas{k});
end

% WRITE OUTPUT VARIABLES
fprintf(fid,'%% OUTPUT\n');
for k = 1:length(variableNames),
   fprintf(fid,'variableValues(%d) = %s;\n',k,variableNames{k});
end 
for k = 1:length(reactionNames),
   fprintf(fid,'reactionValues(%d) = %s;\n',k,reactionNames{k});
end 
fprintf(fid,'return\n');
fprintf(fid,'\n');

% WRITE MODEL FUNCTIONS
fprintf(fid,'%% FUNCTIONS\n');
[filenames,functionFormulas,functionArguments] = SBfunctions(sbm);
for k = 1:length(filenames),
    fprintf(fid,'function [result] = %s(%s)\n',filenames{k},functionArguments{k});
    fprintf(fid,'result = %s;\n',functionFormulas{k});
    fprintf(fid,'return\n');
    fprintf(fid,'\n');
end

% WRITE THE MATLAB FUNCTIONS
fprintf(fid,'%% MATLAB FUNCTIONS\n');
functionsMATLAB = SBfunctionsMATLAB(sbm);
if ~isempty(functionsMATLAB),
    fprintf(fid,'%s',functionsMATLAB);
    fprintf(fid,'\n');
end

% CLOSE FILE
fclose(fid);
return