function [step,dd,delta,newVAL] = brentLineSB(FUN,X,h,p,I,val,dd,delta,jMax)
% Derivative free line Search using quadrtaic interpolation
%
% USAGE:
% ~~~~~~
% [step,dd,delta,newVAL] = brentLineSB(FUN,X,h,p,I,val,dd,delta,jMax)
%
% FUN: Function to minimize in brent
% X: Current guess of the minimizer of FUN
% h: Scale of the distance between X and the minimizer
% p: Search direction
% I: Index of search direction or 0 for extrapolation along a valley
% val: Function value at X ('' if not available)
% dd: Vector holding an appromimative second derivative of the function 
%     phi(lambda)=FUN(X+lambda*U(I)) in dd(I), if not available NaN
% delta(i): amount of decrease of the FUN(x) with the latest steps
%     step(i)*U(:,i)
% jMax: Maximum number of times an attempt is made to halve the trial step
%
% Output Arguments:
% ~~~~~~~~~~~~~~~~~
% step: Approximative minimizing step length
% dd: New calculated or old values of dd
% delta(I): amount of decrease of FUN(x) with step step*U(:,I)
% newVAL: Value of FUN after the taken step step

% Information:
% ~~~~~~~~~~~~
% Systems Biology Toolbox for MATLAB
% Copyright (C) 2007 FCC, mats.jirstrand@fcc.chalmers.se
% Author: Andrea Kirsch
% 
% This program is free software; you can redistribute it and/or modify it
% under the terms of the GNU General Public License as published by the
% Free Software Foundation; either version 2 of the License, or (at your
% option) any later version.
% 
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
% Public License for more details. 
% 
% You should have received a copy of the GNU General Public License along 
% with this program; if not, write to the Free Software Foundation, Inc., 
% 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


global FunEval

% number of variables
ndim = length(X);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TWO INTERPOLATION NODES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nodes and corresponding function values
if I == 0
    % phi(lambda) = FUN(extrapol(lambda)) extrapolation along a valley
    x1 = h(2);                      % = dist(2)
    phi_x0 = p(1);                  % = ZVAL(1) = feval(FUN,Z(:,1))
    phi_x1 = p(3);                  % = ZVAL(3) = feval(FUN,Z(:,3))
else
    % phi(lambda) = FUN(X + lambda*p)
    x1 = h;
    if ~isempty(val)
        phi_x0 = val;
    else
        phi_x0 = feval(FUN,X);
        FunEval = FunEval+1;
    end
    phi_x1 = feval(FUN,X + x1*p);
    FunEval = FunEval+1;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTERPOLATION WITH SECOND DERIVATIVE ESTIMATE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if I > 0 && I <= ndim && ~isnan(dd(I)) && dd(I) ~= 0
    % predict minimum: root of the first derivative of the quadratic
    % interpolating polynomial
    lambda = (1/x1)*((phi_x0 - phi_x1)/dd(I) + 0.5*x1^2);
    
    for j = 1:(jMax+1)
        VAL = feval(FUN,X + lambda*p);
        FunEval = FunEval+1;
        if VAL < phi_x0
            step = lambda;
            delta(I) = phi_x0 - VAL;
            newVAL = VAL;
            return;
        else
            if lambda*x1 <= 0
                lambda = lambda/2;
            else
                break;
            end
        end
    end
    if lambda*x1 <= 0
        step = 0;
        delta(I) = 0;
        newVAL = phi_x0;
        return;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INTERPOLATION WITH A THIRD NODE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% is done if dd(I) not available or I == 0 or when the above interpolation
% was not successful
if phi_x1 <= phi_x0
    x2 = 2*x1;
else
    x2 = -x1;
end
if I == 0
    F2 = feval(@extrapol,h,X,x2);
else
    F2 = X + x2*p;
end
phi_x2 = feval(FUN,F2);
FunEval = FunEval+1;

% calculate second divided difference
div_diff = (((phi_x0 - phi_x1)/(-x1)) - ((phi_x1 - phi_x2)/(x1-x2))) /(-x2);
if I > 0 && I <= ndim
    % approximate second derivative for direction U(:,I)
    dd(I) = 2*div_diff;
end

if div_diff == 0
    % no minimum predicion possible then return step = 0
    step = 0;
    newVAL = phi_x0;
    if I > 0 && I <= ndim
    	delta(I) = 0;
    end
    return;
end

% predict minimum: root of the first derivative of the quadratic 
% interpolating polynomial
lambda = 0.5* (((phi_x0 - phi_x1)/(x1 * div_diff)) + x1);

for j = 1:(jMax+1)
    if I == 0
        estMin = feval(@extrapol,h,X,lambda);
    else
        estMin = X + lambda*p;
    end
    VAL = feval(FUN,estMin);
    FunEval = FunEval+1;
    if VAL < phi_x0
        step = lambda;
        newVAL = VAL;
        if I > 0 && I <= ndim
            delta(I) = phi_x0 - VAL;
        end                
        return;
    else
        lambda = lambda/2;
    end
end

% no success: return step = 0
step = 0;
newVAL = phi_x0;
if I > 0 && I <= ndim
    delta(I) = 0;
end
return;
