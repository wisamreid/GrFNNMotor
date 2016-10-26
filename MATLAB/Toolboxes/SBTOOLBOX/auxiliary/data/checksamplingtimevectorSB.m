function [deltaT] = checksamplingtimevectorSB(deltaTvector,tol)
% checksamplingtimevectorSB
% The idea is to check if the variation between the elements in the passed
% vector does not exceed its mean by more than a relative tolerance of
% "tol" percent. If it eceeds, 0 is returned, otherwise the mean value is
% returned. This functionality is used at sevaral places in the toolbox and
% therefore implemented as an auxiliary script.

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

deltaTvector = unique(deltaTvector);
meanDeltaT = mean(deltaTvector);
toleranceTime = meanDeltaT*tol/100;
toleranceCheck = abs(deltaTvector - meanDeltaT) > toleranceTime;
if sum(toleranceCheck) == 0,
    deltaT = meanDeltaT;
else
    deltaT = 0;
end
return