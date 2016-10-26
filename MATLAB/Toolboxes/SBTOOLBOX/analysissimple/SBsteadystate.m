function [varargout] = SBsteadystate(varargin)
% SBsteadystate: determines a steady-state of an SBmodel or an ODE file 
% description created by SBcreateODEfile. This function should only be used 
% for time-independent systems. Otherwise an error will
% occur. The function is able to deal with models having a singular 
% Jacobian, e.g., due to moity conservations.
%
% USAGE:
% ======
% [steadystate,residual,message] = SBsteadystate(model)         
% [steadystate,residual,message] = SBsteadystate(model,initialCondition)         
% [steadystate,residual,message] = SBsteadystate(model,initialCondition,options)
% [steadystate,residual] = SBsteadystate(model)         
% [steadystate,residual] = SBsteadystate(model,initialCondition)         
% [steadystate,residual] = SBsteadystate(model,initialCondition,options)
% [steadystate] = SBsteadystate(model)         
% [steadystate] = SBsteadystate(model,initialCondition)         
% [steadystate] = SBsteadystate(model,initialCondition,options)
%                                   
% model: SBmodel model or name of ODE file (without .m suffix)
% initialCondition: Array with initial conditions around which to start the
%   search for a steady-state. If not given the initial conditions or given
%   as an empty vector the initial values are taken from the SBmodel.
% options: a vector with options for the solver.
%
%          options = [TolFun,tol,MaxIter,Delta]
%          if options is not given, default options are used.
%          if options are given, then a value of 0 indicates default
%          values.
%          TolFun:  tolerance for max element in function evaluation
%                   (default: 1e-12)
%          tol: Tolerance for the determination of the number of algebraic
%               relationships in the model. If the method fails than this
%               might be due to a wrong rank computation and in this case
%               the tolerance should be increased. If set, the same tolerance
%               setting is used when determining the indices of the
%               dependent variables.
%          MaxIter: maximum number of iterations (default: 1000)
%          Delta:   step length for numerical differentiation to obtain 
%                   jacobian (default: 1e-6)
%
% DEFAULT VALUES:
% ===============
% initialCondition: the initial conditions stored in the model
% options: the default values are used. Choosing 0 for a certain option
%   also leads to the use of the default value for this option.
% TolFun: 1e-11
% tol: standard MATLAB tolerance: s = svd(A); tol = max(size(A))*eps(max(s));
% MaxIter: 1000
% DeltaMax: 1e-6
%
% Output Arguments:
% =================
% steadystate: array with steadystate values
% residual: 2-norm of the derivatives at the found steady-state. residual
%           is optional
% message: message about the steady state not able to find and/or message about 
%          found algebraic relationships / moiety conservations.
%          Message as output argument is optional. If omitted, the message
%          will be displayed in the Matlab window.

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEEDED GLOBAL VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global ODEfctname options message

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INITIALIZATION OF MESSAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
message = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK IF SBMODEL OR FILENAME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp('SBmodel',class(varargin{1})),
    % SBmodel
    sbm = varargin{1};
    % Create temporary ODE file
    [ODEfctname, ODEfilefullpath] = SBcreateTempODEfile(sbm);    
else
    % ODEfctname of ODE file
    ODEfctname = varargin{1};
    % empty sbm variable
    sbm = [];
    % Check if file exists
    if exist(strcat(ODEfctname,'.m')) ~= 2,
        error('ODE file could not be found.');
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS VARIABLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 1,
    initialCondition = feval(ODEfctname);
    options = [0 0 0 0]; % default options
elseif nargin == 2,
    initialCondition = varargin{2};
    options = [0 0 0 0]; % default options
elseif nargin == 3,
    initialCondition = varargin{2};
    options = varargin{3};
else
    error('Wrong number of input arguments');
end
% if empty initial condition => default values
if isempty(initialCondition),
    initialCondition = feval(ODEfctname);
end    
% initialCondition should be a column vector!
initialCondition = initialCondition(:); 
% if options = [] then use default values
if isempty(options),
    options = [0 0 0 0];
end
% check length of options
if length(options) ~= 4,
    error(sprintf('If options are specified you need to give all of them [TolFun,tol,MaxIter,Delta].\nUse 0 values to choose default values for options.')); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET THE STEADY STATE FOR THE MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
steadystate = getSteadyState(initialCondition,sbm);
if ~isempty(steadystate),
    residual = norm(feval(ODEfctname,[],steadystate));
else
    messageText = 'Steady state could not be found. Try different options and/or a different starting guess.';
    message = sprintf('%s\n%s\n',message,messageText);
    residual = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS VARIABLE OUTPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargout == 2,
    varargout{1} = steadystate;
    varargout{2} = residual;
    disp(message);
elseif nargout == 3,
    varargout{1} = steadystate;
    varargout{2} = residual;
    varargout{3} = message;
elseif nargout < 2, 
    varargout{1} = steadystate;
    disp(message);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DELETE FILE IF SBMODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp('SBmodel',class(varargin{1})),
    deleteTempODEfileSB(ODEfilefullpath);
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET THE STEADY STATE OF THE MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [steadystate] = getSteadyState(initialCondition,sbm)
% Define some global variables that are needed in other functions
global depVar factors constantValues ODEfctname options message
TolFun = options(1);
tol = options(2);
MaxIter = options(3);
Delta = options(4);
EXITFLAG = 0;
% Before determining the steady-state it is necessary to check if there are
% some algebraic relationships in the model. If there are not the
% steady-state is determined directly. If there are, the algebraic
% relations are determined and the steady-state is calculated for a reduced
% system.
% check if model given as SBmodel:
if ~isempty(sbm),
    % check if stoichiometric matrix gives full information
    N = SBstoichiometry(sbm,0,1); % no rawFlag,  but silentFlag
    if length(initialCondition) == size(N,1),
        % yes, full information available
        useStoich = 1;
    else
        useStoich = 0;
    end
else
    useStoich = 0;
end
% handle stoich and jac differently
if useStoich == 0,
    % use Jacobian!
    % Get jacobian and determine its rankdeficiency - use default deltaX for
    % numerical differentiation
    J = SBjacobian(ODEfctname,initialCondition);
    if tol == 0,
        d = length(J) - rank(J);        
    else
        d = length(J) - rank(J,tol);
    end
else
    % use stoichiometry
    if tol == 0,
        X = rref([N,eye(size(N,1))]);
    else
        X = rref([N,eye(size(N,1))],tol);
    end        
    Nx = X(:,1:size(N,2));
    d = length(find(sum(abs(Nx)')'==0));
end
if d == 0,
    % Jacobian has full rank => no algebraic relations => "simple"
    % computation of the steady-state
    % The systems does not contain algebraic relations
    [steadystate,FVAL,EXITFLAG] = fsolveSB(@model,initialCondition,[MaxIter,TolFun,Delta]);
else
    % message
    messageText = sprintf('The system has a rank deficiency of %d. This might be due to moiety or\nother conservations in the model, to the excessive use of zero initial\nconditions and/or integrating behavior of the system. If this does\nnot seem to be correct, please check the initial conditions and/or\nset a different value for the tolerance "tol" in the options.',d);
    message = sprintf('%s\n%s\n',message,messageText);
    % Jacobian is singular => assume that algebraic relations are present
    % and solve for the steady-state by taking into account these algebraic
    % relations, by considering the singular output directions
    % Handle stoich and jac differently
    if useStoich == 0,
        message = sprintf('%s\nThe algebraic relationships are determined using the Jacobian.\nThis is numerically not very stable and you need eventually\nto change tolerance "tol" settings.\n',message);
        % Use the singular value decomposition of the Jacobian to determine
        % the algebraic relationships
        [U,S,V] = svd(J);
        U0 = U(:,end-d+1:end);
        if tol == 0,
            factors = rref(U0');
        else
            factors = rref(U0',tol);
        end
    else
        message = sprintf('%s\nThe algebraic relationships are determined using the Stoichiometric matrix.\n',message);
        % Use rref of the Stoichiometric Matrix
       factors = X(end-d+1:end,size(N,2)+1:end);
    end
    constantValues = factors*initialCondition;
    % determine the indices of the dependent variables. this is simple. 
    % just take the index of the first non zero element in each row as 
    % index for dependent a variable.
    depVar = [];
    for k = 1:size(factors,1),
        nonZeroIndices = find(factors(k,:)~=0);
        depVar = [depVar nonZeroIndices(1)];
    end
    % now do the solving for the steady-state
    [ssVarInDep,FVAL,EXITFLAG] = fsolveSB(@modelIndependentVariables,getIndependentVariables(initialCondition),[MaxIter,TolFun,Delta]);
    % Obtain the full state vector from the independent variables
    steadystate = getAllVariables(ssVarInDep);
end
% If fsolve has not converged to a solution then return an empty vector
if EXITFLAG ~= 1,
    steadystate = [];
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REMOVE THE TIME FROM THE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xdot] = model(x)
    global ODEfctname 
    xdot = feval(ODEfctname,[],x);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CALL MODEL BASED ON INDEPENDENT VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [xdotVarInDep] = modelIndependentVariables(varInDep)
    global ODEfctname 
    variables = getAllVariables(varInDep);
    xdot = feval(ODEfctname,[],variables);
    xdotVarInDep = getIndependentVariables(xdot);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET INDEPENDENT VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [varInDep] = getIndependentVariables(variables)
    global depVar
    varInDep = variables(setxor([1:length(variables)],depVar));
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET ALL VARIABLES FROM INDEPENDENT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [variables] = getAllVariables(varInDep)
    global depVar factors constantValues
    n = size(factors,2);
    m = length(depVar);
    varInDepIndex = setxor(1:n,depVar);
    A = factors(1:m,depVar);
    B = constantValues - factors(1:m,varInDepIndex)*varInDep;
    depVarValues = linsolve(A,B);
    variables = zeros(n,1);
    variables(varInDepIndex) = varInDep;
    variables(depVar) = depVarValues;
return      
