function [output] = lipsolSB(A,b,c,lbounds,ubounds,varargin)
% lipsolSB: Linear programming Interior-Point Solver, which solves the
% following type of problem: 
%
%        min           c'*x
%        s.t.         Ax = b
%             lbounds <= x <= ubounds
%
% USAGE:
% ======
% [output] = lipsolSB(A,b,c,lbounds,ubounds)      
% [output] = lipsolSB(A,b,c,lbounds,ubounds,OPTIONS)      
%
% A:        see above
% b:        see above
% c:        see above
% lbounds:  lower bounds for the elements in x. All elements need to be
%           initiailized with reasonable values.
% ubounds:  upper bounds for the elements in x. No upper bound can be
%           specified by choosing 'inf' for an element. In the algorithm
%           this will however be exchanged against the BIG number
%           BIG=1.0000e+032
% OPTIONS: optional argument, allowing the specification of options:
%       OPTIONS.display:    =0 no output to screen
%                           =1 output of iteration information
%                           =2 same as 1 with additional graphical output
%       OPTIONS.maxiter:    maximum number of iterations
%       OPTIONS.cachesize:  cachesize in kB (only change if you know what
%                           you are doing)
%
% DEFAULT VALUES:
% ===============
% OPTIONS.display:      0 (no output to screen)
% OPTIONS.maxiter:      99  
% OPTIONS.cachesize:    1024
%
% Output Arguments:
% =================
% The output is a structure with the following content
% output.xsol:          
% output.results.info	
% output.results.x      
% output.results.y      
% output.results.z      
% output.results.s      
% output.results.w      
%
% Information:
% ============
% Systems Biology Toolbox for MATLAB
% Copyright (C) 2005-2007 Fraunhofer-Chalmers Research Centre
% Author of the toolbox: Henning Schmidt, henning@hschmidt.de
%
% Based on code by Yin Zhang, January, 1995 Department of Mathematics and Statistics
% University of Maryland  Baltimore County modified March, 2003, G. Sweriduk
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
global DISPLAY_LIPSOL
global MAXITER_LIPSOL
global CACHE_SIZE  
global LOOP_LEVEL  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 5,
    OPTIONS = [];
elseif nargin == 6,
    OPTIONS = varargin{1};
else
    error('Incorrect number input of arguments.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set default values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DISPLAY_LIPSOL = 0;
MAXITER_LIPSOL = 99;
CACHE_SIZE = 1024;
LOOP_LEVEL = 8;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Handle OPTIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% display
if isfield(OPTIONS,'display'),
    if ~isempty(OPTIONS.display),
        DISPLAY_LIPSOL = OPTIONS.display;
    end
end
% maxiter
if isfield(OPTIONS,'maxiter'),
    if ~isempty(OPTIONS.maxiter),
        MAXITER_LIPSOL = OPTIONS.maxiter;
    end
end
% maxiter
if isfield(OPTIONS,'cachesize'),
    if ~isempty(OPTIONS.cachesize),
        CACHE_SIZE = OPTIONS.cachesize;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocess ubounds (exchange inf against a big number used in lipsol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ubounds(find(ubounds==inf)) = 1.0000e+032;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Preprocess the problem and check feasibility
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A,b,c,lbounds,ubounds,FEASIBLE] = preprocess(A,b,c,lbounds,ubounds);
if ~FEASIBLE
   error('Infeasibility detected in preprocessing');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Apply scaling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[A,b,c,ubounds] = scaling(A,b,c,ubounds);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine and locate dense columns of matrix A (global variables)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
checkdense(A);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Run the solver
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x,y,z,s,w,info] = miip(A,b,c,ubounds);  % actual solver

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Postprocess result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[xsol, objp] = postprocess(x,lbounds);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generate output argument
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
output = [];
output.xsol = xsol;
output.objp = objp; % c'*xsol
output.converged = info(1);
output.iterations = info(2);
output.trerror = info(3);
output.results.x = x;
output.results.y = x;
output.results.z = x;
output.results.s = x;
output.results.w = x;

