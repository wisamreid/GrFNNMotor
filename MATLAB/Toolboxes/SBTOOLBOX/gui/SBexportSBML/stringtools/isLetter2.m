function [isLetter] = isLetter2(character)
%isLetter2
% tests wether a given character is a Letter
% (it's a simpler and faster implementation for stringhandling)
%
% USAGE:
% ======
%
% [boolean]= isLetter2(character)
%
% character: character to test
%
% boolean: true if character is A, ..., Z, a, ..., z, "_"
%          otherwise false
%
% For more details take a look at the code, please!

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

    % Option 1: letterMatrix = ['a':'z','A':'Z','_'];
    %           if strfind(letterMatrix, character),
    %
    % the latter variant could be usefull if reactionnames can consist of
    % more special characters i.e. ('#', '_', '%', '$', ...)
    % 
    %at the moment this one is the fastest
    if isstrprop(character, 'alpha') || (character == '_'),
        isLetter = 1;
    else
        isLetter = 0;
    end

return