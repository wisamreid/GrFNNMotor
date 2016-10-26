function [X,FVAL] = bfgsSB(FUN,X,varargin)
% BFGS: Quasi-Newton algorithm BFGS method based on section 8.1 in
% "Numerical Optimization", ISBN 0-387-98793-2
%
% USAGE:
% ~~~~~~
% [X,FVAL] = bfgsSB(F,X)
% [X,FVAL] = bfgsSB(F,X,OPTIONS)
%
% FUN: Function to optimize
% X: Starting guess
% OPTIONS: Structure containing options for the algorithm. 
%          OPTIONS.GRAD: Gradient of the function F (as column vector)
%          OPTIONS.TOL: Convergence tolerance
%          OPTIONS.c1: Wolfe condition constant for the line serach to
%               ensure sufficient decrease
%          OPTIONS.c2: Wolfe condition constant for the line search to
%               rule out unacceptable small steps
%          OPTIONS.a1: First trial step length
%          OPTIONS.const: Constant factor to increase step length every
%               iteration of the line search
%          OPTIONS.aMax: Maximum step length
%          OPTIONS.MaxFunEval: Maximum number of function and gradient
%               evaluations
%          OPTIONS.display: Displaying performance: 1 = on, 0 = off
% 
% DEFAULT VALUES:
% ~~~~~~~~~~~~~~~
% OPTIONS.GRAD = @apprGradSB (Approximation of the gradient)
% OPTIONS.TOL  = 1e-6
% OPTIONS.c1   = 1e-4
% OPTIONS.c2   =  0.9
% OPTIONS.a1   =    1
% OPTIONS.const = 1.5
% OPTIONS.aMax  = 500
% OPTIONS.MaxFunEval = 500*numberVariables
% OPTIONS.display = 1
%
% Output Arguments:
% ~~~~~~~~~~~~~~~~~
% X: Found solution
% FVAL: Value of the function F at X

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
% GLOBAL VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global NoLineS NoZoom FunEval GradEval subopt display

% counter for number of lineSearch iterations
NoLineS = 0;
% counter for number of zoom iterations
NoZoom = 0;
% variable to stop the algorithm in case that the line search gives
% suboptimal step lengths several times in a row without effecting a
% progress in the optimization
subopt = 0;
% number of variables
ndim = length(X);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLE OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 2
    OPTIONS = [];
elseif nargin > 2
    OPTIONS = varargin{1};
else
    error('Incorrect number of input arguments.')
end

% default values
GRAD = @apprGradSB;
TOL = 1e-6;
c1 = 1e-4;
c2 = 0.9;
a1 = 1;
const = 1.5;
aMax = 500;
MaxFunEval = 500*ndim;
display = 1;

% read out options
if ~isempty(OPTIONS)
    if isfield(OPTIONS,'GRAD'),
      if ~isempty(OPTIONS.GRAD), GRAD = OPTIONS.GRAD; end
    end
    if isfield(OPTIONS,'TOL'),
        if ~isempty(OPTIONS.TOL), TOL = OPTIONS.TOL; end
    end
    if isfield(OPTIONS,'c1'),
        if ~isempty(OPTIONS.c1), c1 = OPTIONS.c1; end
    end
    if isfield(OPTIONS,'c2'),
        if ~isempty(OPTIONS.c2), c2 = OPTIONS.c2; end
    end
    if isfield(OPTIONS,'a1'),
        if ~isempty(OPTIONS.a1), a1 = OPTIONS.a1; end
    end
    if isfield(OPTIONS,'const'),
        if ~isempty(OPTIONS.const), const = OPTIONS.const; end
    end
    if isfield(OPTIONS,'aMax'),
        if ~isempty(OPTIONS.aMax), aMax = OPTIONS.aMax; end
    end
    if isfield(OPTIONS,'MaxFunEval'),
        if ~isempty(OPTIONS.MaxFunEval), MaxFunEval = OPTIONS.MaxFunEval; end
    end
    if isfield(OPTIONS,'display'),
        if ~isempty(OPTIONS.display), display = OPTIONS.display; end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZATION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% avoiding dimension conflicts
R = size(X,1);
if R == 1
    X = X';
end
x_old = X;
% initial inverse Hessian approximation
H = eye(ndim);

% test variable to memorise in which mode the program operates:
%   apprGradSBient = 1, gradient is calculated approximately
%   apprGradSBient = 0, gradient is given explicitly as input
apprGradSBient = 1;
if ~isequal(GRAD,@apprGradSB)
    apprGradSBient = 0;
end

FVAL = feval(FUN,X);
if apprGradSBient
    grad = feval(GRAD,x_old,FUN,TOL);
else
    grad = feval(GRAD,x_old);
end

% counter for number of BFGS iterations
BFGSiter = 0;
% counter for number of function and gradient evaluations
FunEval = 1;
GradEval = 1;

% displaying performance
if display 
    disp('No BFGS Iter   No LineS Iter    No Zoom Iter   No Fun Eval    No Grad Eval   Min Function   Step Length   Algo Step');
    if apprGradSBient
        disp(sprintf('%5u %14u %16u %14u %14u %21.6g %29s', 0, 0, 0, FunEval+2*ndim*GradEval, 0, FVAL, 'Initial Guess'));
    else
        disp(sprintf('%5u %14u %16u %14u %14u %21.6g %29s', 0, 0, 0, FunEval, GradEval, FVAL, 'Initial Guess'));
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% test variable to know if H has to be scaled
scaling = 1;
FVAL_old = inf;

while norm(grad) > TOL
    % search direction
    p = -H*grad;
    % compute step length in lineSearch
    alpha = bfgslineSB(FUN,GRAD,x_old,FVAL,p,TOL,apprGradSBient,c1,c2,a1,const,aMax);
    if alpha == inf
        break;
    end
    % new guess of minimizer
    x_new = x_old + alpha*p;
    if apprGradSBient
        grad_new = feval(GRAD,x_new,FUN,TOL);
    else
        grad_new = feval(GRAD,x_new);
    end    
    GradEval = GradEval+1;
    s = x_new - x_old;
    y = grad_new - grad;
    % making the size of H similar to that of the inverse Hessian, only
    % during the first interation or when the Hessian couldn't be computed
    % before and was set to the identity matrix instead
    x_old = x_new;
    grad = grad_new;
    FVAL = feval(FUN,x_old);
    FunEval = FunEval+1;
    if scaling
        % exception handling:
        if y'*y == 0
            % neither the BFGS update nor setting H to a scaled identity
            % matrix can be done to find a reasonable serach direction
            if display 
                if apprGradSBient
                    disp(sprintf('%5u %14u %16u %14u %14u %21.6g', BFGSiter, NoLineS, NoZoom, FunEval+2*ndim*GradEval, 0, FVAL));
                else
                    disp(sprintf('%5u %14u %16u %14u %14u %21.6g', BFGSiter, NoLineS, NoZoom, FunEval, GradEval, FVAL));
                end
            end
            disp('Divide by zero. H cannot be computed.');
            break; 
        end
        H = ((y'*s)/(y'*y))*H;
        scaling = 0;
    end
    if (y'*s) == 0
        % then H cannot be updated by the BFGS formula
        % instead set H = eye(ndim) and scale it after the next step length
        % has been computed like in the beginning
        if norm(grad_new) <= TOL
            % then H does not have to be updated anymore
            break;
        elseif FVAL >= (FVAL_old - 0.0001)
            % trial with identity matrix did not effect a satisfying
            % progress
            disp('Divide by zero. H cannot be computed and setting to eye instead did not effect a satisfying progress either.');
            break;
        else
            H = eye(ndim);
            scaling = 1;
            % memorize the current function value to check later if
            % setting H to the identity matrix was reasonable
            FVAL_old = FVAL;
            if display
                if apprGradSBient
                    disp(sprintf('%5u %14u %16u %14u %14u %21.6g  %33s', BFGSiter, NoLineS, NoZoom, FunEval+2*ndim*GradEval, 0, FVAL, 'Resetting H to Eye'));
                else
                    disp(sprintf('%5u %14u %16u %14u %14u %21.6g  %33s', BFGSiter, NoLineS, NoZoom, FunEval, GradEval, FVAL, 'Resetting H to Eye'));
                end
            end
        end
    else
        % BFGS update for inverse approximative Hessian matrix
        rho = 1/(y'*s);
        H = (eye(ndim) - rho*s*y')*H*(eye(ndim) - rho*y*s') + rho*s*s';
    end
    % counting and displaying performance
    BFGSiter = BFGSiter+1;
    if display && scaling == 0
        if apprGradSBient
            disp(sprintf('%5u %14u %16u %14u %14u %21.6g  %26s', BFGSiter, NoLineS, NoZoom, FunEval+2*ndim*GradEval, 0, FVAL, 'BFGS Update'));
        else
            disp(sprintf('%5u %14u %16u %14u %14u %21.6g  %26s', BFGSiter, NoLineS, NoZoom, FunEval, GradEval, FVAL, 'BFGS Update'));
        end
    end
    if apprGradSBient
        if FunEval + 2*ndim*GradEval > MaxFunEval
            disp('Exceeded maximum number of function evaluations.');
            break;
        end
    else
        if FunEval + GradEval > MaxFunEval
            disp('Exceeded maximum number of function and gradient evaluations.');
            break;
        end
    end
end
% the output equals the latest guess of the minimizer x_new. In case the
% while loop is never run it equals the starting guess.
X = x_old;
return
