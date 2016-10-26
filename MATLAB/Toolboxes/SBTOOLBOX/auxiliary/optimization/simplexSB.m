function [X,FVAL,EXITFLAG] = simplexSB(FUN,X,varargin)
% simplexSB: Downhill Simplex Method in Multidimensions. Algorithm based on
% section 10.4 in "Numerical Recipes in C", ISBN 0-521-43108-5.
% The algorithm has been modified in order to be able to take into account
% simple box constraints on parameter values.
% 
% USAGE:
% ======
% [X,FVAL,EXITFLAG] = simplexSB(FUN,X,OPTIONS)
%
% FUN: Function to optimize
% X: Starting Guess
% OPTIONS: structure containing options for the algorithm:
%        OPTIONS.maxfunevals: Maximum number of function evaluations
%        OPTIONS.maxiter: Maximum number of iterations
%        OPTIONS.maxtime: Maximum time (in minutes) for optimization
%        OPTIONS.tolfun: Tolerance for difference between best and worst
%                        function evaluation in simplex
%        OPTIONS.tolx: Tolerance for max difference between the coordinates of
%                      the vertices.
%        OPTIONS.highbounds: vector containing upper bounds for parameters
%           Instead of a vector highbounds can also be a scalar > 1. In the
%           latter case the highbounds are determined by multiplying this
%           scalar to the initial parameter guess X. This latter case works
%           only for positive X(i) ... if negative present then the user
%           needs to specify a whole lowbounds vector. 
%        OPTIONS.lowbounds: vector containing lower bounds for parameters.
%           Instead of a vector lowbounds can also be a scalar < 1. In the
%           latter case the lowbounds are determined by multiplying this
%           scalar to the initial parameter guess X. This latter case works
%           only for positive X(i) ... if negative present then the user
%           needs to specify a whole highbounds vector. 
%        OPTIONS.outputFunction: string with output function name. If
%               not given or if empty then no output function will be used.
%               This output function can be used to display data, or to 
%               save optimization information to a file. The function needs
%               to be in the MATLAB path. Its calling syntax is:
%                   'outputFunction'(bestparameters,bestfunctionvalue,currentsimplex)
%        OPTIONS.silent: =0: output of info, =1: no output
%
% DEFAULT VALUES:
% ===============
% OPTIONS.maxfunevals:    200*numberVariables
% OPTIONS.maxiter:        200*numberVariables
% OPTIONS.maxtime:        120
% OPTIONS.tolx:           1e-10
% OPTIONS.tolfun:         1e-10
% OPTIONS.lowbounds:      0.1  => lowbounds = 0.1*X 
% OPTIONS.highbounds:     10  => highbounds = 10*X 
% OPTIONS.outputFunction: no output function ('')
% OPTIONS.silent:         0 (no output of info)
%
% Output Arguments:
% =================
% X: Found solution
% FVAL: Value of the function FUN at X
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
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,
% USA.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLOBAL VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global ndim nfunk Xguess lowbounds highbounds

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 2,
    OPTIONS = [];
elseif nargin == 3,
    OPTIONS = varargin{1};
else
    error('Incorrect number input of arguments.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
EXITFLAG = 1;
ndim = length(X);
maxfunevals = 200*ndim;
maxiter = 200*ndim;
maxtime = 120;
tolx = 1e-10;
tolfun = 1e-10;
outputFunction = '';
lowbounds = 0.1*X;
highbounds = 10*X;
silent = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Handle OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% silent
if isfield(OPTIONS,'silent'),
    if ~isempty(OPTIONS.silent),
        silent = OPTIONS.silent;
    end
end
% tolfun
if isfield(OPTIONS,'tolfun'),
    if ~isempty(OPTIONS.tolfun),
        tolfun = OPTIONS.tolfun;
    end
end
% tolx
if isfield(OPTIONS,'tolx'),
    if ~isempty(OPTIONS.tolx),
        tolx = OPTIONS.tolx;
    end
end
% maxiter
if isfield(OPTIONS,'maxtime'),
    if ~isempty(OPTIONS.maxtime),
        maxtime = OPTIONS.maxtime;
    end
end
% maxiter
if isfield(OPTIONS,'maxiter'),
    if ~isempty(OPTIONS.maxiter),
        maxiter = OPTIONS.maxiter;
    end
end
% maxfunevals
if isfield(OPTIONS,'maxfunevals'),
    if ~isempty(OPTIONS.maxfunevals),
        maxfunevals = OPTIONS.maxfunevals;
    end
end
% lowbounds
if isfield(OPTIONS,'lowbounds'),
    if ~isempty(OPTIONS.lowbounds),
        if length(OPTIONS.lowbounds) == ndim,
            lowbounds = OPTIONS.lowbounds;
        else
            if length(OPTIONS.lowbounds) == 1,
                if ~isempty(find(X<0)),
                    error('Please specify a vector for the lower bounds using OPTIONS.lowbounds.');
                end
                lowbounds = X*OPTIONS.lowbounds;
            else
                error('The OPTIONS.lowbounds setting is not correct.');
            end
        end
    end
else
    if ~isempty(find(X<0)),
        error('Please specify a vector for the lower bounds using OPTIONS.lowbounds.');
    end
end
% highbounds
if isfield(OPTIONS,'highbounds'),
    if ~isempty(OPTIONS.highbounds),
        if length(OPTIONS.highbounds) == ndim,
            highbounds = OPTIONS.highbounds;
        else
            if length(OPTIONS.highbounds) == 1,
               if ~isempty(find(X<0)),
                    error('Please specify a vector for the higher bounds using OPTIONS.highbounds.');
               end
                highbounds = X*OPTIONS.highbounds;
            else
                error('The OPTIONS.highbounds setting is not correct.');
            end
        end
    end
else
    if ~isempty(find(X<0)),
        error('Please specify a vector for the higher bounds using OPTIONS.highbounds.');
    end
end
lowbounds = lowbounds(:)';
highbounds = highbounds(:)';
% outputFunction
outputFunction = '';
if isfield(OPTIONS,'outputFunction'),
    if ~isempty(OPTIONS.outputFunction),
        outputFunction = OPTIONS.outputFunction;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK BOUND VIOLATION FOR INITIAL GUESS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
indexXhi = find(X(:) > highbounds(:));
if ~isempty(indexXhi),
    error('Initial guess does violate high parameter bounds.');
end
indexXlo = find(X(:) < lowbounds(:));
if ~isempty(indexXlo),
    error('Initial guess does violate low parameter bounds.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZE EMPTY SIMPLEX DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
p = zeros(ndim+1,ndim);     % vertice vectors in rows
y = zeros(ndim+1,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FIRST POINT OF INITIAL SIMPLEX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xguess = X(:)'; 
p(1,:) = Xguess;    
y(1) = feval(FUN,Xguess);

if ~silent,
    disp(' Nr Iter    Nr Fun Eval         Min function          Algorithm Step');
    disp(sprintf(' %5.0f        %5.0f           %12.6g            %s', 0, 1, y(1), 'initial guess'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REMAINING POINTS OF INITIAL SIMPLEX
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% construct vertices by modifying one element each
% relative changes in case that elements are non-zero,
% absolute changes in case that elements are zero
relativeDelta = 0.05; 
absoluteDelta = 0.00025;
for k = 1:ndim
    Xmodify = Xguess;
    if Xmodify(k) == 0
        % absolute change
        if highbounds(k) > absoluteDelta,
            Xmodify(k) = absoluteDelta;
        else 
            Xmodify(k) = -absoluteDelta;
        end
    else
        % relative change
        Xmodify(k) = (1 + relativeDelta)*Xmodify(k);
    end
    p(k+1,:) = Xmodify;
    y(k+1) = feval(FUN,Xmodify);
end
algostep = 'initial simplex';
nriterations = 1;
nfunk = ndim+1;

% if output function given then run output function to plot
% intermediate result
if length(outputFunction) ~= 0,
    feval(outputFunction,p(1,:), y(1),p);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN ALGORITHM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reorder y and p so that the first row corresponds to the 
% lowest function value
tic % start timer
while(1), 
    % do sorting instead of determining the indices of the best, worst,
    % next worst
    help = sortrows([y,p],1);
    y = help(:,1);
    p = help(:,2:end);
    if ~silent,
        disp(sprintf(' %5.0f        %5.0f           %12.6g            %s', nriterations, nfunk, y(1), algostep));
    end
    % if output function given then run output function to plot
    % intermediate result
    if length(outputFunction) ~= 0,
        feval(outputFunction,p(1,:), y(1),p);
    end
    % end the optimization if the difference between best and worst
    % function evaluation in simplex is smaller than tolfun and the
    % max difference between the coordinates of the verttices is less than
    % tolx
    if abs(y(end)-y(1)) < tolfun && max(max(abs(p(2:ndim+1)-p(1:ndim)))) < tolx,
        break;
    end
    % check number of function evaluations
    if nfunk >= maxfunevals,
        EXITFLAG = 0;
        if ~silent,
            disp('Exceeded maximum number of function evaluations.');
        end
        break;
    end
    % check number of iterations
    if nriterations >= maxiter,
        EXITFLAG = 0;
        if ~silent,
            disp('Exceeded maximum number of iterations.');
        end
        break;
    end
    if toc/60 > maxtime,
        EXITFLAG = 0;
        if ~silent,
            disp('Exceeded maximum time.');
        end
        break;
    end
    % Begin a new iteration. First extrapolate by a factor -1 through the face of the simplex
    % across from the high point, i.e., reflect the simplex from the high point.
    [ytry,ptry] = amotry(FUN, p, y, -1);
    % check the result 
    if ytry <= y(1),
        % Gives a result better than the best point, so try an additional 
        % extrapolation by a factor 2.
        [ytryexp,ptryexp] = amotry(FUN, p, y, -2);
        if ytryexp < ytry,
            p(end,:) = ptryexp;
            y(end) = ytryexp;
            algostep = 'extrapolation';
        else
            p(end,:) = ptry;
            y(end) = ytry;
            algostep = 'reflection';
        end
    elseif ytry >= y(ndim),
        % The reflected point is worse than the second-highest, so look 
        % for an intermediate lower point, i.e., do a one-dimensional
        % contraction. 
        [ytrycontr,ptrycontr] = amotry(FUN, p, y, -0.5);
        if ytrycontr < y(end),
            p(end,:) = ptrycontr;
            y(end) = ytrycontr;
            algostep = 'one dimensional contraction';
        else
            % Can’t seem to get rid of that high point. Better contract
            % around the lowest (best) point.
            x = ones(ndim,ndim)*diag(p(1,:));
            p(2:end,:) = 0.5*(p(2:end,:)+x);
            for k=2:ndim,
                y(k) = feval(FUN,p(k,:));
            end
            algostep = 'contraction around best point';
        end
    else
        % if ytry better than second-highest point then use this point
        p(end,:) = ptry;
        y(end) = ytry;
        algostep = 'reflection';
    end        
    nriterations = nriterations + 1;
end
% do the sorting a last time
help = sortrows([y,p],1);
y = help(:,1);
p = help(:,2:end);
X = p(1,:);
FVAL = y(1);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AMOTRY FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ytry,ptry] = amotry(FUN, p, y, fac)
% Extrapolates by a factor fac through the face of the simplex across from 
% the high point, tries it, and replaces the high point if the new point is 
% better.
global ndim nfunk lowbounds highbounds
psum = sum(p(1:ndim,:))/ndim;
ptry = psum*(1-fac) + p(end,:)*fac;

% Deal with low and high parameter bounds
indexXhi = find(ptry > highbounds);
indexXlo = find(ptry < lowbounds);
ptry(indexXhi) = highbounds(indexXhi);
ptry(indexXlo) = lowbounds(indexXlo);

% Evaluate the function at the trial point.
ytry = feval(FUN,ptry);
nfunk = nfunk + 1;
return

