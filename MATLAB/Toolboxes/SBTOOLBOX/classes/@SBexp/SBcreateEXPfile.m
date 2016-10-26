function SBcreateEXPfile(varargin)
% SBcreateEXPfile: creates a *.exp file with the experiments text description
%
% USAGE:
% ======
% [] = SBcreateEXPfile(exp)         
% [] = SBcreateEXPfile(exp,filename)
%
% exp: SBexp object to convert to a textfile description
% filename: filename for the created textfile 
%
% DEFAULT VALUES:
% ===============
% filename: constructed from the experiments name

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
if nargin == 1,
    exp = varargin{1};
    % if no filename provided then use the name of the SBexp object as filename
    % remove unwanted characters first
    functionName = regexprep(exp.name,'\W','');
    filename = strcat(functionName,'.exp');
elseif nargin == 2,
    exp = varargin{1};
    filename = strcat(varargin{2},'.exp');
else
    error('Incorrect number of input arguments.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONVERT MODEL TO TEXT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[expTextStructure] = convertExpToTextSB(exp);
[completeText] = setPartsToCompleteTextExpSB(expTextStructure);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITE THE FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(filename,'w');
fwrite(fid,completeText);
fclose(fid);
return
