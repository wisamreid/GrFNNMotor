function [stepLength] = bfgslineSB(FUN,GRAD,X,VAL,p,TOL,apprGradient,c1,c2,a1,const,aMax)
% lineSearch: Line search Algorithm for the Wolfe conditions based on
% section 3.4 in "Numerical Optimization", ISBN 0-387-98793-2
%
% USAGE:
% ~~~~~~
% [stepLength] = bfgslineSB(FUN,GRAD,X,p,c1,c2,a1,const,aMax)
% 
% FUN: function to optimize in BFGS
% GRAD: Gradient of the function FUN
% X: Current guess of the minimizer of FUN in BFGS
% VAL: function value at X
% p: Search direction
% TOL: Error tolerance for the function apprGradSB
% apprGradSBient: mode in which the program operates, =1 if GRAD=@apprGradSB
% c1: Wolfe condition constant to ensure sufficient decrease in the
%     objective function (e.g. c1=1e-4)
% c2: Wolfe condition constant to rule out unacceptable small steps
%     (e.g. c2=0.9 when the search direction is chosen by a Newton or
%     quasi-Newton method)
% a1: First trial step length
% const: Constant factor to increase step length every iteration
% aMax: Maximum step length
%
% Output Arguments:
% ~~~~~~~~~~~~~~~~~
% stepLength: Found step length holding the Wolfe conditions or aMax

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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZATION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global NoLineS NoZoom FunEval GradEval subopt display

ndim = length(X);
a_old = 0;
a_new = a1;
phi0 = VAL;
phiVal_old = phi0;
if apprGradSBient
    phiGrad0 = feval(GRAD,X,FUN,TOL);
else
    phiGrad0 = feval(GRAD,X);
end
GradEval = GradEval+1;
i = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while(1)
    NoLineS = NoLineS+1;
    phiVal_new = feval(FUN,X + a_new*p);
    FunEval = FunEval+1;
    if phiVal_new > phi0 + c1*a_new*phiGrad0'*p || (i>1 && phiVal_new >= phiVal_old)
        stepLength = zoom(a_old,phiVal_old,a_new,FUN,GRAD,X,p,TOL,apprGradSBient,phi0,phiGrad0,c1,c2);
        if stepLength == inf
            break;
        end        
        % displaying performance
        if display 
            if apprGradSBient
                disp(sprintf('%20u %16u %14u %14u %35.3g   %s', NoLineS, NoZoom, FunEval+2*ndim*GradEval, 0, stepLength, 'LineSearch'));
            else
                disp(sprintf('%20u %16u %14u %14u %35.3g   %s', NoLineS, NoZoom, FunEval, GradEval, stepLength, 'LineSearch'));
            end
        end
        break;
    end
    if apprGradSBient
        phiGrad_new = feval(GRAD,X + a_new*p,FUN,TOL);
    else
        phiGrad_new = feval(GRAD,X + a_new*p);
    end    
    GradEval = GradEval+1;
    if abs(phiGrad_new'*p) <= -c2*phiGrad0'*p
        stepLength = a_new;
        subopt = 0;
        % displaying performance
        if display 
            if apprGradSBient
                disp(sprintf('%20u %16u %14u %14u %35.3g   %s', NoLineS, NoZoom, FunEval+2*ndim*GradEval, 0, stepLength, 'LineSearch'));
            else
                disp(sprintf('%20u %16u %14u %14u %35.3g   %s', NoLineS, NoZoom, FunEval, GradEval, stepLength, 'LineSearch'));
            end
        end
        break;
    end
    if phiGrad_new'*p >= 0
        stepLength = zoom(a_new,phiVal_new,a_old,FUN,GRAD,X,p,TOL,apprGradSBient,phi0,phiGrad0,c1,c2);
        if stepLength == inf
            break;
        end
        % displaying performance
        if display 
            if apprGradSBient
                disp(sprintf('%20u %16u %14u %14u %35.3g   %s', NoLineS, NoZoom, FunEval+2*ndim*GradEval, 0, stepLength, 'LineSearch'));
            else
                disp(sprintf('%20u %16u %14u %14u %35.3g   %s', NoLineS, NoZoom, FunEval, GradEval, stepLength, 'LineSearch'));
            end
        end
        break;
    end
    % update trial step length, here as constant multiple of the latest
    % one, but could also be done by means of interpolation procedures
    a_help = const*a_new;
    if a_help > aMax
        if display 
            disp('Exceeded maximum step length in lineSearch. Therefore aMax was chosen as step length.');
        end
        stepLength = aMax;
        subopt = 0;
        % displaying performance
        if display 
            if apprGradSBient
                disp(sprintf('%20u %16u %14u %14u %35.3g   %s', NoLineS, NoZoom, FunEval+2*ndim*GradEval, 0, stepLength, 'LineSearch'));
            else
                disp(sprintf('%20u %16u %14u %14u %35.3g   %s', NoLineS, NoZoom, FunEval, GradEval, stepLength, 'LineSearch'));
            end
        end
        break
    end
    % memorize the function value at a_new which is a_old in the next
    % iteration
    phiVal_old = phiVal_new;
    a_old = a_new;
    a_new = a_help;
    i = i+1;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ZOOM FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [alpha] = zoom(a_lo,VALa_lo,a_hi,FUN,GRAD,X,p,TOL,apprGradSBient,phi0,phiGrad0,c1,c2)
% Finds acceptable step length alpha in an interval that brackets step
% lengths which hold the Wolfe conditions.

global NoZoom FunEval GradEval subopt display

ndim = length(X);
j = 0;
while(1)
    NoZoom = NoZoom+1;
    % interpolation to find trial step length a_j between a_lo and a_hi,
    % here by means of bisection, but could also be done using quadratic or
    % cubic interpolation
    a_j = (a_lo+a_hi)/2;
    phiVal_j = feval(FUN,X + a_j*p);
    FunEval = FunEval+1;
    phiVal_lo = VALa_lo;
    if phiVal_j > phi0 + c1*a_j*phiGrad0'*p || phiVal_j >= phiVal_lo
        a_hi = a_j;
        if phiVal_lo == phiVal_j
            j = j+1;
        end
        if j > 10
            if phiVal_j == phi0
                % in case this exceptional choice of the step length is no
                % good stategie
                subopt = subopt + 1;
                if subopt > 10
                    disp('Using suboptimal step lengths is no promising stategy.');
                    alpha = inf;
                    break;
                end
            end
            if display 
                disp('Exceeded 10 trial step lengths without attaining a lower function value. Latest suboptimal trial was accepted.');
            end
            alpha = a_j;
            % displaying performance
            if display 
                if apprGradSBient
                    disp(sprintf('%37u %14u %14u %35.3g   %s', NoZoom, FunEval+2*ndim*GradEval, 0, alpha, 'Zoom'));
                else
                    disp(sprintf('%37u %14u %14u %35.3g   %s', NoZoom, FunEval, GradEval, alpha, 'Zoom'));
                end
            end
            break;
        end
    else
        if apprGradSBient
            phiGradVal = feval(GRAD,X + a_j*p,FUN,TOL);
        else
            phiGradVal = feval(GRAD,X + a_j*p);
        end        
        GradEval = GradEval+1;
        if abs(phiGradVal'*p) <= -c2*phiGrad0'*p
            alpha = a_j;
            subopt = 0;
            % displaying performance
            if display 
                if apprGradSBient
                    disp(sprintf('%37u %14u %14u %35.3g   %s', NoZoom, FunEval+2*ndim*GradEval, 0, alpha, 'Zoom'));
                else
                    disp(sprintf('%37u %14u %14u %35.3g   %s', NoZoom, FunEval, GradEval, alpha, 'Zoom'));
                end
            end
            break;
        end
        if phiGradVal'*p*(a_hi-a_lo) >= 0
            a_hi = a_lo;
        end
        a_lo = a_j;
        j = 0;
    end
end
return
