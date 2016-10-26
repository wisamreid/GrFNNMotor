function [formula] = formula2vecSBPD(formula)
% formula2vecSBPD: Converts a formula given as a string into a formula that
% can be evaluated with variables given as vectors. (Adding '.' in front of
% '*', '/', '^'

% Information:
% ============
% SBPD Package - Systems Biology Parameter Determination Package
% based on SBaddon code (c) 2006 Henning Schmidt, FCC
% 
% Contact: Henning Schmidt, henning@sbtoolbox.org
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

formula = regexprep(formula,'\*','.*');
formula = regexprep(formula,'\/','./');
formula = regexprep(formula,'\^','.^');

return