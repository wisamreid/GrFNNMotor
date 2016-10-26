function [output] = SBsensdatastat(varargin)
% SBsensdatastat: Function allowing to generate steady-state data that 
% subsequently can be used for different kinds of parametric sensitivity 
% analyses. 
%
% The data is obtained by determining the steady-states of the nominal 
% and perturbed systems. In each calculation of a perturbed system only 
% a single parameter is perturbed.
% As initial guess for the steady-state the initial conditions stored in the
% model are used. As nominal parameter values the parameter values stored
% in the model are used.
% In case the considered steady-state is unstable at least a good guess of
% the steady-state needs to be given in the model.
%  
% USAGE:
% ======
% [output] = SBsensdatastat(model)
% [output] = SBsensdatastat(model,parameters)
% [output] = SBsensdatastat(model,parameters,pertSize,absRel)
% [output] = SBsensdatastat(model,parameters,pertSize,absRel,options)
%
% model: SBmodel (not useable with ODE model file)
% parameters: a cell-array with parameter names, determining the parameters
%       that are to be considered in the analysis.
% pertSize: either a scalar value, determining the perturbation to be
%       applied to each of the parameters, or a vector of same length as
%       the number of considered parameters, allowing to choose a different 
%       perturbation for each parameter.
% absRel: either a scalar (0 or 1), or a vector of the same length as
%       pertSize. 0 means that the corresponding element in pertSize
%       determines an absolute perturbatiom, 1 means it determines a
%       relative perturbation in percent.
% options: options are passed as they are to the SBsteadystate function.
%          options = [TolFun,tol,MaxIter,Delta] more information in the
%          help text of the SBsteadystate function.
%
% DEFAULT VALUES:
% ===============
%   The only required input argument is 'model'.
%   If not specified otherwise, following default values are used:
%
%   parameters = all parameters in the model
%   pertSize = 1 (percent)
%   absRel = 1 (relative perturbation)
%   options = [] (default options of the SBsteadystate function)
%
% OUTPUT:
% =======
% The output argument 'output' is a MATLAB struct element with the
% following structure:
%
%   output.model            model as obtained as input argument
%   output.states           cell-array with names of states as elements
%   output.xssnom           vector containing the steady-state values of
%                           the states of the nominal system 
%   output.xsspert          cell-array, where each entry correponds to 
%                           a vector containing the steady-state of the
%                           states of the perturbed systems
%   output.reactions        cell-array with names of reactions as elements
%   output.rssnom           vector containing the steady-state values of
%                           the reaction rates of the nominal system 
%   output.rsspert          cell-array, where each entry correponds to 
%                           a vector containing the steady-state of the
%                           reaction rates of the perturbed systems
%   output.parameters       same as the input argument 'parameters' -
%                           either its default or the passed values
%   output.nomvalues        nominal values of parameters in above field
%   output.pertSize         same as the input argument 'pertSize' -
%                           either its default or the passed values
%   output.absRel           same as the input argument 'absRel' -
%                           either its default or the passed values
%
% The reason of using a struct element as output is that the output of 
% SBsensdatastat usually needs to be processed by other functions, and the
% calling of these functions is considerably simplified by passing all
% important data by one element only.

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
% CHECK SBMODEL 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sbm = varargin{1};
if ~strcmp('SBmodel',class(sbm)),
    error('This function can only be used with SBmodels, not with ODE file models!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET NAME FOR ODE FILE AND MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% path to systems temporary directory
TEMPDIR = tempdir;  
% add path to temporary directory
addpath(TEMPDIR);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% HANDLING VARIABLE INPUT ARGUMENTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 1,
    parameters = SBparameters(sbm);
    pertSize = ones(1,length(parameters));
    absRel = ones(1,length(parameters));
    states = SBstates(sbm);
    options = [];
elseif nargin == 2,
    parameters = varargin{2};
    if strcmp('char',class(parameters)),
        % adjusting for the case where just one parameter is given as a
        % string
        parameters = {parameters};
    end
    pertSize = ones(1,length(parameters));
    absRel = ones(1,length(parameters));
    states = SBstates(sbm);
    options = [];
elseif nargin == 4,
    parameters = varargin{2};
    pertSize = varargin{3};
    absRel = varargin{4};
    states = SBstates(sbm);
    options = [];
elseif nargin == 5,
    parameters = varargin{2};
    pertSize = varargin{3};
    absRel = varargin{4};
    states = SBstates(sbm);
    options = varargin{5};
else
    error('Incorrect number of input arguments.');
end
% check if parameters defined by a cell-array 
% if not then convert
if strcmp(class(parameters),'char'),
    parameters = {parameters};
end
% adjust scalar elements in pertSize and absRel to number of parameters
if length(pertSize) == 1,
    pertSize = pertSize*ones(1,length(parameters));
end
if length(absRel) == 1,
    absRel = absRel*ones(1,length(parameters));
end
% check length of parameters, pertSize and absRel
if length(parameters) ~= length(pertSize) || length(parameters) ~= length(absRel),
    error(sprintf('Check the number of elements in ''parameters'', ''pertSize'', and ''absRel''\ninput arguments. They should have the same number of elements!'));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SETTING OTHER VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% determining the nominal values for the parameters
[allParameters, allParameterValues] = SBparameters(sbm);
nomValues = [];
errorText = '';
for k1 = 1:length(parameters),
    parameterFound = 0;
    for k2 = 1:length(allParameters),
        if strcmp(parameters{k1},allParameters{k2}),
            nomValues = [nomValues, allParameterValues(k2)];
            parameterFound = 1;
            break;
        end
    end
    if parameterFound == 0,
        errorText = sprintf('%sParameter ''%s'', given in input arguments could not be found in the model.\n',errorText,parameters{k1});
    end
end
if ~isempty(errorText),
    errorText = sprintf('%sThe available parameters in the model are:\n',errorText);
    for k = 1:length(allParameters),
        errorText = sprintf('%s\n%s',errorText,allParameters{k});
    end    
    error(errorText);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INFO MESSAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Collecting data for sensitivity analysis by steady-state calculations.');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine steady-state for nominal model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
disp('Steady-state computation for nominal system');
[xnom,residual,message] = SBsteadystate(sbm,SBinitialconditions(sbm),options);
if isempty(xnom),
    error('No steady-state could be found. Please try another starting guess.');
end
[dummy1,dummy2,dummy3,rnom] = SBreactions(sbm,xnom);
rnom(find(abs(rnom)<eps)) = 0;
time_elapsed = toc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine steady-states for perturbed models
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Steady-state computation for perturbed systems');
time_estimated = (time_elapsed*length(parameters))/60;
disp(sprintf('Estimated time for calculations: %5.1f minutes.\n',time_estimated));
xpert = {};
rpert = {};
for k = 1:length(parameters),
    % determine the absolute perturbation for current parameter
    if absRel(k) == 0 || nomValues(k) == 0,
        % absolute perturbation
        if nomValues(k) == 0,
            pertParamValue = nomValues(k) + pertSize(k)/100;
        else
            pertParamValue = nomValues(k) + pertSize(k);
        end
        disp(sprintf('%d) Absolute perturbation (%+g) of parameter ''%s''. Nominal: %g, perturbed: %g',k,pertSize(k),parameters{k},nomValues(k),pertParamValue));
        if nomValues(k) == 0,
            disp(sprintf('\tNominal value of parameter ''%s'' is zero. Using absolute perturbation instead of relative.',parameters{k}));
            absRel(k) = 0;
            pertSize(k) = pertSize(k)/100;
        end
    else
        % relative perturbation
        pertParamValue = nomValues(k) * (1 + pertSize(k)/100);
        disp(sprintf('%d) Relative perturbation (%+g%%) of parameter ''%s''. Nominal: %g, perturbed: %g',k,pertSize(k),parameters{k},nomValues(k),pertParamValue));
    end
    % sbm is the original model - determine a perturbed model modelpert -
    % and create an ODE file for it. rehash the path
    sbmpert = SBparameters(sbm,parameters{k},pertParamValue);
    % get steady-state for perturbed system
    [xpertk,residual,message] = SBsteadystate(sbmpert,SBinitialconditions(sbmpert),options);
    xpert{k} = real(xpertk);
    [dummy1,dummy2,dummy3,rpertk] = SBreactions(sbmpert,xpert{k});
    rpertk(find(abs(rpertk)<eps)) = 0;
    rpert{k} = rpertk;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCT OUTPUT STRUCTURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
output = [];
output.model = sbm;
output.states = states;
output.xssnom = xnom;
output.xsspert = xpert;   
output.reactions = SBreactions(sbm);     
output.rssnom = rnom;          
output.rsspert = rpert;         
output.parameters = parameters;  
output.nomvalues = nomValues;
output.pertSize = pertSize;      
output.absRel = absRel;          
return
