function SBcreateODEfile(varargin)
% SBcreateODEfile: creates an m-file that can be simulated by 
% using standard MATLAB integrators (ODE45, ODE23s, etc.)
%
% USAGE:
% ======
% [] = SBcreateODEfile(sbm)         
% [] = SBcreateODEfile(sbm,filename)
% [] = SBcreateODEfile(sbm,filename,dataFileFlag)
% [] = SBcreateODEfile(sbm,filename,dataFileFlag,eventFlag)
%
% sbm: SBmodel to convert to an ODE file description
% filename: filename for the created ODE file (without extension)
% dataFileFlag: =1: creates an additional file allowing to determine
%                   variable and reaction rate values for given state values 
%                   and time vector.
%               =0: doe not create the additional file
% eventFlag: =1: creates 2 additional files for the handling of events in the model. 
%            =0: does not create these two files.
%            THE EVENT FILES ARE ONLY CREATED IN CASE THAT EVENTS ARE
%            PRESENT IN THE MODEL
%            
% DEFAULT VALUES:
% ===============
% filename: constructed from the models name
% dataFileFlag: 0
% eventFlag: 0

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
% PROCESS VARIABLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dataFileFlag = 0;
eventFlag = 0;
if nargin == 1,
    sbm = varargin{1};
    % if no filename provided then use the name of the SBmodel as filename
    % remove white spaces from the name.
    functionName = strrep(sbm.name,' ','');
    functionName = strrep(functionName,'-','');
    functionName = strrep(functionName,'+','');
    functionName = strrep(functionName,'*','');
    functionName = strrep(functionName,'/','');
    filename = strcat(functionName,'.m');
elseif nargin == 2,
    sbm = varargin{1};
    % extract filename from input argument number 2 - needed since 
    % SBsimulate gives a whole path as filename
    [PATHSTR,functionName,EXT,VERSN] = fileparts(varargin{2});
    if length(PATHSTR) == 0,
        filename = strcat(functionName,'.m');
    else
        filename = strcat(PATHSTR,'/',functionName,'.m');
    end
elseif nargin == 3,
    sbm = varargin{1};
    % extract filename from input argument number 2 - needed since 
    % SBsimulate gives a whole path as filename
    [PATHSTR,functionName,EXT,VERSN] = fileparts(varargin{2});
    if length(PATHSTR) == 0,
        filename = strcat(functionName,'.m');
    else
        filename = strcat(PATHSTR,'/',functionName,'.m');
    end
    dataFileFlag = varargin{3};
elseif nargin == 4,
    sbm = varargin{1};
    % extract filename from input argument number 2 - needed since 
    % SBsimulate gives a whole path as filename
    [PATHSTR,functionName,EXT,VERSN] = fileparts(varargin{2});
    if length(PATHSTR) == 0,
        filename = strcat(functionName,'.m');
    else
        filename = strcat(PATHSTR,'/',functionName,'.m');
    end
    dataFileFlag = varargin{3};
    eventFlag = varargin{4};
else
    error('Incorrect number of input arguments.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK IF MODEL CONTAINS STATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(sbm.states) == 0,
    error('The model does not contain any states');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLE THE EVENT FLAG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if the eventflag is set to a non-zero value AND at least on event is
% present in the model an event function file is created that has to be 
% used when calling the integrator. additionally a function is created 
% that determines the new values to which the states have to be set. 
if eventFlag == 1 && length(sbm.events) > 0,
    filenameEvent = strrep(filename,functionName,strcat('event_',functionName));
    filenameEventAssignment = strrep(filename,functionName,strcat('event_assignment_',functionName));
    SBcreateEventFunction(sbm,filenameEvent);
    SBcreateEventAssignmentFunction(sbm,filenameEventAssignment);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLE THE DATAFILE FLAG
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if dataFileFlag == 1,
    filenameDataFile = strrep(filename,functionName,strcat('datafile_',functionName));
    SBcreateSimulationDataFunction(sbm,filenameDataFile);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITE WARNING MESSAGE IN CASE EVENT FLAG IS NOT SET AND EVENTS ARE PRESENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if eventFlag == 0 && length(sbm.events) > 0,
    if length(sbm.events) == 1,
        disp(sprintf('Warning: 1 event is present in the model.\nFor this analysis it is not taken into account.'));
    else
        disp(sprintf('Warning: %d events are present in the model.\nFor this analysis they are not taken into account.',length(sbm.events)));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPEN FILE FOR WRITING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(filename,'w');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITE FUNCTION DEFINITION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'function [output] = %s(varargin)\n',functionName);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITE THE HEADER
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'%% %s\n',sbm.name);
fprintf(fid,'%% Generated: %s\n',datestr(now));
fprintf(fid,'%% \n');
fprintf(fid,'%% [output] = %s() => output = initial conditions in column vector\n',functionName);
fprintf(fid,'%% [output] = %s(''states'') => output = species names in cell-array\n',functionName);
fprintf(fid,'%% [output] = %s(''parameters'') => output = parameter names in cell-array\n',functionName);
fprintf(fid,'%% [output] = %s(''parametervalues'') => output = parameter values in column vector\n',functionName);
fprintf(fid,'%% [output] = %s(time,statevector) => output = time derivatives in column vector\n',functionName);
fprintf(fid,'%% \n');
fprintf(fid,'%% State names and ordering:\n');
fprintf(fid,'%% \n');
for k = 1:length(sbm.states),
    fprintf(fid,'%% statevector(%d): %s\n',k,sbm.states(k).name);
end
fprintf(fid,'%% \n');
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLE VARARGINS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'%% HANDLE VARIABLE INPUT ARGUMENTS\n');
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'if nargin == 0,\n');
    % Return initial conditions for state variables
    fprintf(fid,'\t%% Return initial conditions of the state variables\n');
    if length(sbm.states) > 0,
        fprintf(fid,'\toutput = [');
    else
        fprintf(fid,'\toutput = [];\n');
    end
    count = 0;
    for k = 1:length(sbm.states)-1,
        fprintf(fid,'%s, ',num2str(sbm.states(k).initialCondition,20));
        count = count + 1;
        if count == 10, fprintf(fid,'...\n\t\t'); count = 0; end
    end
    if length(sbm.states) > 0,
        fprintf(fid,'%f];\n',sbm.states(length(sbm.states)).initialCondition);
    end
    fprintf(fid,'\toutput = output(:);\n');
    fprintf(fid,'\treturn\n');
fprintf(fid,'elseif nargin == 1,\n');
	fprintf(fid,'\tif strcmp(varargin{1},''states''),\n');
        fprintf(fid,'\t\t%% Return state names in cell-array\n');
        if length(sbm.states) > 0,
            fprintf(fid,'\t\toutput = {');
        else
            fprintf(fid,'\t\toutput = {};\n');
        end 
        count = 0;
        for k = 1:length(sbm.states)-1,
            fprintf(fid,'''%s'', ',sbm.states(k).name);
            count = count + 1;
            if count == 10, fprintf(fid,'...\n\t\t\t'); count = 0; end
        end
        if length(sbm.states) > 0,
            fprintf(fid,'''%s''};\n',sbm.states(length(sbm.states)).name);
        end
    fprintf(fid,'\telseif strcmp(varargin{1},''parameters''),\n');
        fprintf(fid,'\t\t%% Return parameter names in cell-array\n');
        if length(sbm.parameters) > 0,
            fprintf(fid,'\t\toutput = {');
        else
            fprintf(fid,'\t\toutput = {};\n');
        end 
        count = 0;
        for k = 1:length(sbm.parameters)-1,
            fprintf(fid,'''%s'', ',sbm.parameters(k).name);
            count = count + 1;
            if count == 10, fprintf(fid,'...\n\t\t\t'); count = 0; end
        end
        if length(sbm.parameters) > 0,
            fprintf(fid,'''%s''};\n',sbm.parameters(length(sbm.parameters)).name);
        end
    fprintf(fid,'\telseif strcmp(varargin{1},''parametervalues''),\n');
        fprintf(fid,'\t\t%% Return parameter values in column vector\n');
        if length(sbm.parameters) > 0,
            fprintf(fid,'\t\toutput = [');
        else
            fprintf(fid,'\t\toutput = [];\n');
        end 
        count = 0;
        for k = 1:length(sbm.parameters)-1,
            fprintf(fid,'%s, ',num2str(sbm.parameters(k).value, 20));
            count = count + 1;
            if count == 10, fprintf(fid,'...\n\t\t\t'); count = 0; end
        end
        if length(sbm.parameters) > 0,
            fprintf(fid,'%s];\n',num2str(sbm.parameters(length(sbm.parameters)).value,20));
        end
    fprintf(fid,'\telse\n');
        fprintf(fid,'\t\terror(''Wrong input arguments! Please read the help text to the ODE file.'');\n');
    fprintf(fid,'\tend\n');
    fprintf(fid,'\toutput = output(:);\n');
    fprintf(fid,'\treturn\n');
fprintf(fid,'elseif nargin == 2,\n');
    fprintf(fid,'\ttime = varargin{1};\n');
    fprintf(fid,'\tstatevector = varargin{2};\n');
fprintf(fid,'elseif nargin == 3,\n');
    fprintf(fid,'%% For the reason of initialization and finishing, NARGIN is sometimes 3.\n');
    fprintf(fid,'%% The third input argument id of no interest for this toolbox!\n');
    fprintf(fid,'\ttime = varargin{1};\n');
    fprintf(fid,'\tstatevector = varargin{2};\n');
fprintf(fid,'else\n');
    fprintf(fid,'\terror(''Wrong input arguments! Please read the help text to the ODE file.'');\n');
fprintf(fid,'end\n');
fprintf(fid,'\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE THE STATES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'%% STATES\n');
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
for k = 1:length(sbm.states)
    fprintf(fid,'%s = statevector(%d);\n',sbm.states(k).name,k);
end
fprintf(fid,'\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE THE PARAMETES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(sbm.parameters),
    fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf(fid,'%% PARAMETERS\n');
    fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    for k = 1:length(sbm.parameters),
        fprintf(fid,'%s = %s;\n',sbm.parameters(k).name,num2str(sbm.parameters(k).value,20));
    end
    fprintf(fid,'\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE THE VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(sbm.variables),
    fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf(fid,'%% VARIABLES\n');
    fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    for k = 1:length(sbm.variables),
        fprintf(fid,'%s = %s;\n',sbm.variables(k).name,sbm.variables(k).formula);
    end
    fprintf(fid,'\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE THE REACTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(sbm.reactions),
    fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf(fid,'%% REACTION KINETICS \n');
    fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    for k = 1:length(sbm.reactions),
        fprintf(fid,'%s = %s;\n',sbm.reactions(k).name,sbm.reactions(k).formula);
    end
    fprintf(fid,'\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITE THE ODES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'%% DIFFERENTIAL EQUATIONS\n');
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
for k = 1:length(sbm.states),   
    fprintf(fid,'%s_dot = %s;\n',sbm.states(k).name,sbm.states(k).ODE);
end
fprintf(fid,'\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITE RETURN VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
fprintf(fid,'%% RETURN VALUES\n');
fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
for k = 1:length(sbm.states),   
    fprintf(fid,'output(%d) = %s_dot;\n',k,sbm.states(k).name);
end
fprintf(fid,'%% return a column vector \n');
% Return a column vector
fprintf(fid,'output = output(:);\n');
fprintf(fid,'return\n');
fprintf(fid,'\n');
fprintf(fid,'\n');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITE THE FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(sbm.functions),
    fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf(fid,'%% FUNCTIONS\n');
    fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    for k = 1:length(sbm.functions),
        fprintf(fid,'function [result] = %s(%s)\n',sbm.functions(k).name,sbm.functions(k).arguments);
        fprintf(fid,'result = %s;\n',sbm.functions(k).formula);
        fprintf(fid,'return\n');
        fprintf(fid,'\n');
    end
    fprintf(fid,'\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITE THE MATLAB FUNCTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(sbm.functionsMATLAB),
    fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf(fid,'%% MATLAB FUNCTIONS\n');
    fprintf(fid,'%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n');
    fprintf(fid,'%s',sbm.functionsMATLAB);
    fprintf(fid,'\n');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CLOSE THE FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
status = fclose(fid);
% Return
return
