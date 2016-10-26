function [result] = piecewiseSB(varargin)
% piecewiseSB: This function implements support for the SBML / MATHML
% piecewise operator.
% 
% USAGE:
% ======
% [result] = piecewiseSB(resultiftrue1,decision1,resultiftrue2,decision2,...,resultiftruen,decisionn)   
% [result] = piecewiseSB(resultiftrue1,decision1,resultiftrue2,decision2,...,resultiftruen,decisionn,defaultresult)    
%
% decision1,...,decisionn: logical argument, e.g. returned from a comparison
% result1,...,resultn: returnvalue in case the corresponding decision is
%   evaluated to be true
% defaultresult: if none of the decisions are true this defaultresult is returned.
%
% Output Arguments:
% =================
% result: the result corresponding to the decision that is true. if no
%   decision is true and no defaultresult i given an error will occurr.

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
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
% USA.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DETERMINE THE OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
result = [];
% check if odd or even number of input arguments
oddnumber = mod(nargin,2);
for k = 1:2:nargin-oddnumber,
    if varargin{k+1},
        result = varargin{k};
        break;
    end
end
if isempty(result),
    if oddnumber,
        result = varargin{nargin};
    else
        error('A piecewise statement is wrongly defined - missing (but needed) default value.');
    end
end
return
   