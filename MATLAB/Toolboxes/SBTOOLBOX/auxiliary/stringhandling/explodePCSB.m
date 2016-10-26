function [elements] = explodePCSB(text,varargin)
% explodePCSB: This function does not(!!!) lead to an explosion of your
% Personal Computer. It is an auxiliary function allowing to decompose a
% string expression into its comma separated elements. Commas within
% parentheses expressions are not considered. This is a useful function to
% obtain the arguments of a function call, where some arguments might have
% expressions involving commas but have parentheses around them.
% Alternatively, a separator character, other then a comma, can be
% specified by the user.
%
% USAGE:
% ======
% [elements] = explodePCSB(text)
% [elements] = explodePCSB(text,separatorCharacter)
%
% text: text to decompose into its comma separated elements
% separatorCharacter: one character that should be used for the explosion
% of the string.
%
% DEFAULT VALUES:
% ===============
% separatorCharacter: ',' (comma)
%
% Output Arguments:
% =================
% elements: cell-array containing string elements, previously separated by
% the separator character.

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
% HANDLE VARIABLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 1,
    separatorCharacter = ',';
elseif nargin == 2,
    separatorCharacter = varargin{1};
else
    error('Incorrect number of input arguments.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DO THE EXPLOSION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elements = {};
openParenthesis = 0;
lastIndex = 1;
elementIndex = 1;
for k2 = 1:length(text),
    if text(k2) == '(',
        openParenthesis = openParenthesis + 1;
    end
    if text(k2) == ')',
        openParenthesis = openParenthesis - 1;
    end
    if text(k2) == separatorCharacter && openParenthesis == 0,
        elements{elementIndex} = strtrim(text(lastIndex:k2-1));
        elementIndex = elementIndex + 1;
        lastIndex = k2+1;
    end
end
elements{elementIndex} = strtrim(text(lastIndex:end));
return

