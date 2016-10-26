function [wellFormed] = isBracketstringWellformed(string)
%isBracketStringWellformed
% checks wether the order of opened and closed braces as well as
% their count is valid
%
% USAGE:
% ======
%
% [wellFormed] = isBracketstringWellformed(string)
%
% string: i.e. a formula where correct bracketing is to test
%
% wellFormed: true or false

% Information:
% ============
% Systems Biology Toolbox for MATLAB
% Copyright (C) 2005-2007 Fraunhofer-Chalmers Research Centre
% Author of the toolbox: Henning Schmidt, henning@hschmidt.de
% Programmed by Gunnar Drews, Student at University of
% Rostock Department of Computerscience, gunnar.drews@uni-rostock.de
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

    wellFormed = 0;
    openedBrackets={'{','[','('};
    closedBrackets={'}',']',')'};
    stack='';
    index=1;
    for k = 1 : length(string),
        if ~isempty(strmatch(string(k), openedBrackets)),
            stack(index) = string(k);
            % stack
            index = index + 1;
        elseif ~isempty(strmatch(string(k), closedBrackets)),
            if isempty(stack),
                return;
            elseif (int8(string(k))-1 == int8(stack(index-1))) || (int8(string(k))-2 == int8(stack(index-1))),
                index = index - 1;
                stack = stack(1:index-1);
                % stack
            else
                return;
            end
        end
    end
    if isempty(stack),
        wellFormed = 1;
    end
return