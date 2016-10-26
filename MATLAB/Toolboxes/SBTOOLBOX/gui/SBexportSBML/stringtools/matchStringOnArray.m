function [rowPos] = matchStringOnArray(patternString, matchingArray)
%matchStringOnArray
% A string match function which gives back the row number of the
% first occurence of the pattern string within the matching array.
% (it's a simpler and faster implementation for stringhandling)
%
% USAGE:
% ======
% [rowPos] = matchStringOnArray(patternString, matchingArray)
%
% matchingArray: an array of String
% patternString: a string to test for occurance in the array
%
% rowPos: 0 if string doesn't occur in array
%         >0 (row number) otherwise -> can be used as boolean value

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

rowPos = 0;
searchResult = strcmp(patternString, matchingArray);
for n1 = 1 : length(searchResult),
    if searchResult(n1),
        rowPos = n1;
        break;
    end
end
return