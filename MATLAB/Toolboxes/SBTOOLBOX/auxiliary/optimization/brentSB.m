function [X,FVAL] = brentSB(FUN,X,H,varargin)
% Brent: Principal Axis Method for Minimization without Derivatives
% based on chapter 7 in "Algorithms for Minimization without Derivatives",
% ISBN-10: 0-13-022335-2
%
% USAGE:
% ~~~~~~
% [X,FVAL] = brentSB(FUN,X,H)
% [X,FVAL] = brentSB(FUN,X,H,OPTIONS)
%
% FUN: Function to optimize
% X: Starting guess
% H: Estimate of the distance to the minimum of FUN
% OPTIONS: Structure containing options for the algorithm.
%          OPTIONS.TOL: Positive absolute tolerance
%          OPTIONS.jMax: Maximum number of times an attempt is made to
%               halve the trial step in the line search
%          OPTIONS.pos: Small positive number for the diagonal matrix D
%               used instead of a calculated negative divided difference
%          OPTIONS.stop: Number of consecutive iterations for which the
%               stopping criterion is to be satisfied
%          OPTIONS.MaxFunEval: Maximum number of function evaluations
%          OPTIONS.illc: Scalar value 1 or 0. 1 = Problem known as
%               ill-conditioned, 0 = not
%          OPTIONS.display: Displaying performance: 1 = on, 0 = off
%
% DEFAULT VALUES:
% ~~~~~~~~~~~~~~~
% OPTIONS.TOL  = 1e-5
% OPTIONS.jMax = 5
% OPTIONS.pos  = 0.001;
% OPTIONS.stop = 2
% OPTIONS.MaxFunEval = 1000*numberVariables
% OPTIONS.illc = 0 (problem not ill-conditioned)
% OPTIONS.display = 1 (on)
%
% Output Arguments:
% ~~~~~~~~~~~~~~~~~
% X: Found solution
% FVAL: Value of the function FUN at X

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
global FunEval

% number of variables
ndim = length(X);
FVAL = feval(FUN,X);
FunEval = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLE OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 3
    OPTIONS = [];
elseif nargin == 4
    OPTIONS = varargin{1};
else
    error('Incorrect number of input arguments.')
end

% default values
TOL = 1e-5;
jMax = 5;
pos = 0.001;
stop = 2;
MaxFunEval = 1000*ndim;
illc = 0;
display = 1;

% read out options
if ~isempty(OPTIONS)
    if isfield(OPTIONS,'TOL'),
      if ~isempty(OPTIONS.TOL), TOL = OPTIONS.TOL; end
    end
    if isfield(OPTIONS,'jMax'),
      if ~isempty(OPTIONS.jMax), jMax = OPTIONS.jMax; end
    end
    if isfield(OPTIONS,'pos'),
      if ~isempty(OPTIONS.pos), pos = OPTIONS.pos; end
    end
    if isfield(OPTIONS,'stop'),
      if ~isempty(OPTIONS.stop), stop = OPTIONS.stop; end
    end
    if isfield(OPTIONS,'MaxFunEval'),
      if ~isempty(OPTIONS.MaxFunEval), MaxFunEval = OPTIONS.MaxFunEval; end
    end
    if isfield(OPTIONS,'illc'),
      if ~isempty(OPTIONS.illc), illc = OPTIONS.illc; end
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
x0 = X;
h = H;

% matrix of initial search directions
U = eye(ndim);
S = NaN*eye(ndim);
% vector containing the an approximate second derivative, at position dd(i)
% the one corresponding to search direction U(:,i), if not available NaN
dd = NaN*ones(ndim,1);
% vector to store the in the line search computed devided differences
D = NaN*ones(ndim,1);
% vector to store the latest step lengths at a time, beta(i) will hold the
% taken step length in direction U(:,i)
beta = NaN*ones(ndim,1);                % how to initialize? book says: if
                                        % beta(i)=0 we use a result of a
                                        % previous iteration, what before?
% amount of decrease of FUN(x) with the latest steps beta(i)*U(:,i)
delta = zeros(ndim,1);
% vector to store the values to be maximized in order to find the best
% discarding direction
disc = -inf*ones(ndim,1);

% matrix to store 3 consecutive minimum approximations in case an
% extrapolbrentSBation along a valley is to be done
Z = NaN*ones(ndim,3);
Z(:,1) = X;
% vector to store the correspondent function values
ZVAL = NaN*ones(1,3);
ZVAL(1) = FVAL;

% several counter
BrentIter = 0;
NoLineS = 0;
NoSVD = 0;
% number of iterations from the start or after an SVD 
iter = 0;
% number of consecutive iterations satisfying the stopping criterion 
STOP = 0;
RESET = 0;
res = 0;

% variable to trigger a random step trial before a new iteration
RS = 0;

% displaying performance
if display
    disp('No Iter   No LineS Iter   No Fun Eval   Min Function   Step Length   Algorithm Step');
    disp(sprintf('%5u %14u %13u %17.6g %29s', BrentIter, NoLineS, FunEval, FVAL, 'Initial Guess'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_new = x0;
while(1)
    iter = iter + 1;
    % not in the first iteration and after resetting U
    if iter ~= 1
        % if the problem is ill-conditioned or if recent linear searches 
        % have failed to improve the current approximation to the minimum
        % try random step to get off resolution valley
        if RS || illc
            % step of length 10*RE in a random direction where RE is about
            % the rounding error in the computation of FUN
            if ~isnan(S(1,1)) && sqrt(2*norm(FVAL)*eps/(max(max(S)))^(-2)) ~= 0
                RE = sqrt(2*norm(FVAL)*eps/(max(max(S)))^(-2));
            else
                RE = TOL;        % how to choose better      ??
            end
            RND = rand(ndim,1);
            % scale RND so that norm(U*RND)=10*RE
            RND = (10*RE/norm(U*RND))*RND;
            x_new = x0 + U*RND;
            if display
                disp(sprintf('%80s','Random Step'));
            end
        else
            RND = zeros(ndim,1);
        end
        for i = 1:ndim
            % search successive in all ndim directions U(:,i)
            if i == 1
                if RS
                    val = '';
                    RS = 0;
                else
                    val = ZVAL(1);                    
                end                
            elseif i == 2 && artStep
                % in case beta(1) was set to h
                val = '';
            else
                val = ValNew;
            end
            [stepLength,dd,delta,ValNew] = brentLineSB(FUN,x_new,h,U(:,i),i,val,dd,delta,2);
            NoLineS = NoLineS+1;
            % ensure that beta(1) ~= 0 so the directions can never become 0
            if i == 1 && stepLength == 0
                stepLength = h;
                artStep = 1;
            else
                artStep = 0;
            end
            if stepLength ~= 0
                % if stepLength == 0 use previous result to find dicarding
                % direction (otherwise division by zero)
                beta(i) = stepLength;
            end
            x_new = x_new + stepLength*U(:,i);
        end
        
        % choosing direction to discard
        for j = 1:(ndim-iter+1)
              disc(j) = (beta(j) + RND(j))*sqrt(delta(j))/norm(beta(j));                            
        end
        % search for the index j which maximizes disc(j)
        [nn,DISC] = max(disc((ndim-iter+1)));
        % overwrite the discarding direction
        U(:,DISC) = U(:,(ndim-iter+1));
        dd(DISC) = dd((ndim-iter+1));
        % store the new direction at the first position that will not be
        % changed in the next iteration (ndim-iter+1)
        U(:,(ndim-iter+1)) = x_new - x0;
        dd((ndim-iter+1)) = NaN;
    end
    
    % search along the new found conjugate direction
    [stepLength,dd,delta,newVAL] = brentLineSB(FUN,x0,h,U(:,(ndim-iter+1)),(ndim-iter+1),ZVAL(1),dd,delta,jMax);
    NoLineS = NoLineS+1;
    x_new = x0 + stepLength*U(:,(ndim-iter+1));
    FVAL = newVAL;
    if stepLength == 0;
        RS = 1;
    else
        % store the new minimum estimation and the correspondent function
        % values for later extrapolbrentSBation along valley, only if it is better
        % (stepLength~=0) 
        for j = 3:-1:2
            Z(:,j) = Z(:,j-1);
            ZVAL(j) = ZVAL(j-1);
        end
        Z(:,1) = x_new;
        ZVAL(1) = FVAL;
    end
    % count and display performance
    BrentIter = BrentIter+1;
    if display
        disp(sprintf('%5u %14u %13u %17.6g %13.3g   %s', BrentIter, NoLineS, FunEval, FVAL, stepLength, 'Brent Iteration'));
    end
    if FunEval > MaxFunEval
        disp('Exceeded maximum number of function evaluations.');
        X = x_new;
        break;
    end
    % new estimate of the distance to the minimum
    if ~isnan(ZVAL(3)) && norm(Z(:,3) - Z(:,1)) < h
        h = norm(Z(:,3) - Z(:,1));
    end    

    if iter == ndim
        if NoSVD >= 2 && ~isnan(ZVAL(3))
            % extrapolbrentSBation along a valley
            dist(1) = norm(Z(:,2) - Z(:,1));
            dist(2) = norm(Z(:,1) - Z(:,3));
            if dist(1) ~= 0 && dist(1) ~= 0
                % otherwise division by zero
                [stepLength,dd,delta,newVAL] = brentLineSB(FUN,Z,dist,ZVAL,0,'',dd,delta,2);
                x_new = extrapolbrentSB(dist,Z,stepLength);
                FVAL = newVAL;
                % store the new found minimum estimation and the 
                % correspondent function values
                if stepLength ~= 0
                    for j = 3:-1:2
                        Z(:,j) = Z(:,j-1);
                        ZVAL(j) = ZVAL(j-1);
                    end
                    Z(:,1) = x_new;
                    ZVAL(1) = FVAL;
                end
                if display
                    disp(sprintf('%52.6g %13.3g %22s', FVAL, stepLength, 'Valley extrapolbrentSBation'));
                end
            end
        end
        
        % resetting U after every ndim iterations to avoid linear
        % dependence of the the search directions so that the search for
        % the minimum cannot become restricted to a subspace
        for k = 1:ndim
            D(k) = dd(k)/2;         % = div_diff(k)
            if D(k) <= 0
                D(k) = pos;    
            end
            if D(k) == inf || isnan(D(k))
                % svd cannot be calculated, try resetting U to eye
                res = 1;
            end
        end
        if res
            RESET = RESET+1;
        else
            RESET = 0;
            res = 0;
        end
        if RESET >= 1 && RESET <= 10
            U = eye(ndim);
            S = NaN*eye(ndim);
            dd = NaN*ones(ndim,1);
            iter = 0;
            if display
                disp(sprintf('%87s', 'Resetting U to Eye'));
            end
        elseif RESET >= 11
            disp('Singular value decomposition cannot be calculated anymore.')
            X = x_new;
            return;
        else
            V = U*((diag(D))^(-0.5));
            % new set U of search directions by singular value decomposition
            [U,S,nn] = svd(V);
            NoSVD = NoSVD+1;
            for j=1:ndim
                % new approximate second derivatives
                dd(j) = 2* (S(j,j))^(-2);
            end
            iter = 0;
            if display
                disp(sprintf('%72s', 'SVD'));
            end
        end
    end
    
    % stopping criterion, if satisfied for the desired number 'stop' of
    % consecutive iterations the algorithm terminates with the
    % approximation x_new
    if 2*norm(x0 - x_new) <= TOL    % 2*norm(x0 - x_new) <= eps^(0.5)*norm(x_new) + TOL
        STOP = STOP+1;
        if STOP >= stop
            X = x_new;
            if display
                disp(sprintf('%5u %14u %13u %17.6g %27s', BrentIter, NoLineS, FunEval, FVAL, 'Termination'));
            end
            break;
        end
    else
        STOP = 0;        
    end
      
    % memorizing the current best approximation of the minimizer
    x0 = x_new;    
end

FVAL = feval(FUN,X);
return
