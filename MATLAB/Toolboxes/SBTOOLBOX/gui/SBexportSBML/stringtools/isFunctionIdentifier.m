function [isFunctionIdentifier] = isFunctionIdentifier(identifier, varargin)
%isFunctionIdentifier
% checks wether identifier is a reserved function identifier
%
% USAGE:
% ======
% [boolean] = isFunctionIdentifier(identifier)
% identifier: test token
%
% boolean: true if identifier is a reserved MATLAB function identifier
%          false otherwise
%
% [boolean] = isFunctionIdentifier(identifier, SBmodel) 
% identifier: test token
% SBmodel: SBmodel 
%
% boolean: true if identifier is a reserved MATLAB function identifier or
%          if identifier is already used in the supplied SBmodel
%          false otherwise
%
% For a list of reserved function identifiers take a look at the code,
% please!

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

modelGiven = false;
if nargin == 2,
    % get all function names from given model
    sbmodel = varargin{1};
    modelGiven = true;
end

    isFunctionIdentifier = false;
    % Corresponding MATLAB expressions
    FunctionExpressions = {'acos', 'asin', 'atan', 'ceil', 'log',...
        'log10', 'power', 'acos', 'acosh', 'acot', 'acoth', 'acsc',...
        'acsch', 'asec', 'asech', 'asin', 'asinh', 'atan', 'atanh',...
        'exp', 'logbase', 'piecewise', 'and', 'or', 'gt', 'lt', 'eq',...
        'ge', 'le', 'mod', 'piecewiseSB', 'andSB', 'orSB'};
    for index = 1 : length(FunctionExpressions),
        if strcmp(identifier, FunctionExpressions{index}),
            isFunctionIdentifier = true;
            return;
        end
    end
    
    % if SBmodel given, test function names there, too
    if modelGiven,
        functionNames = getSBmodelFunctionNames(sbmodel);
        for index = 1 : length(functionNames),
            if strcmp(identifier, functionNames{index}),
                isFunctionIdentifier = true;
                return;
            end
        end
    end
return