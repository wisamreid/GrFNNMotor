function [X,FVAL,EXITFLAG] = fsolveSB(FUN,X,varargin)
% fsolveSB: attempts to solve equations of the form FUN(X)=0    
% where FUN and X may be vectors. Newton iteration.
% 
% USAGE:
% ======
% [X,FVAL,EXITFLAG] = fsolveSB(FUN,X)
% [X,FVAL,EXITFLAG] = fsolveSB(FUN,X,OPTIONS)
%
% FUN: Sunction to minimize/solve
% X: Starting Guess
% OPTIONS: a vector with options for the solver.
%          OPTIONS = [MaxIter,TolFun,Delta]
%          if options is not given, default options are used.
%          if options are given, then a value of 0 indicates default
%          values.
%          MaxIter: maximum number of iterations
%          TolFun:  tolerance for max element in function evaluation
%          Delta:   step length for numerical differentiation to 
%                   obtain jacobian
%
% DEFAULT VALUES:
% ===============
% MaxIter: 1000
% TolFun: 1e-11
% Delta: 1e-6
%
% Output Arguments:
% =================
% X: Found solution
% FVAL: Value of the equations FUN at X
% EXITFLAG: 1=success, 0=not found

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

if nargin == 3,
    OPTIONS = varargin{1};
    if isempty(OPTIONS),
        OPTIONS = [0 0 0];
    else
        if length(OPTIONS) ~= 3,
            error('Incorrect number of elements in options vector.');
        end 
    end
else
    OPTIONS = [0 0 0];
end

if OPTIONS(1) == 0, MaxIter = 1000; else MaxIter = OPTIONS(1); end
if OPTIONS(2) == 0, TolFun = 1e-11; else TolFun = OPTIONS(2); end
if OPTIONS(3) == 0, Delta = 1e-8; else Delta = OPTIONS(3); end


EXITFLAG = 0;
for k = 1:MaxIter,
    FVAL = feval(FUN,X);
    Jacobian = getJacobian(FUN,X,Delta);
    X = X - 0.5*pinv(Jacobian)*FVAL;
    if norm(FVAL) < TolFun,
        EXITFLAG = 1;
        break
    end
end

if EXITFLAG == 0,
    disp('SBsteadystate/fsolveSB: Exceeded maximum number of iterations - check options');
    disp(sprintf('Last residual: %g',max(abs(FVAL))));
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET JACOBIAN AT GIVEN STATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Jacobian] = getJacobian(FUN,X,DELTA)
% size of system
n = length(X);          
% initialize jacobian variable
Jacobian = zeros(n,n);  
% determine the Jacobian by numerical differentiation
for k = 1:n,                
    Xup = X;
    Xup(k) = Xup(k)+DELTA;
    Jacobian(:,k) = (FUN(Xup)'-FUN(X)')/DELTA;
end
return








