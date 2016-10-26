function SBcreateTEXTBCfile(varargin)
% SBcreateTEXTBCfile: creates a *.txtbc file with the models text
% description in a biochemically oriented format
%
% USAGE:
% ======
% [] = SBcreateTEXTBCfile(sbm)         
% [] = SBcreateTEXTBCfile(sbm,filename)
%
% sbm: SBmodel to convert to a textfile description
% filename: filename for the created textfile (.txtbc)
%
% DEFAULT VALUES:
% ===============
% filename: constructed from the models name

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
    sbm = varargin{1};
    % if no filename provided then use the name of the SBmodel as filename
    % remove unwanted characters first
    functionName = regexprep(sbm.name,'\W','');
    filename = strcat(functionName,'.txtbc');
elseif nargin == 2,
    sbm = varargin{1};
    filename = strcat(varargin{2},'.txtbc');
else
    error('Incorrect number of input arguments.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONVERT MODEL TO TEXT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[modelTextBCStructure] = convertModelToTextBCSB(sbm);
[completeTextBC] = setPartsToCompleteTextBCSB(modelTextBCStructure);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% WRITE THE FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fid = fopen(filename,'w');
fwrite(fid,completeTextBC);
fclose(fid);
return
