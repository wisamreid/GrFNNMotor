function [output] = reducereac(output,indexkeep,varargin)
% reducereac: reduces the equation by keeping only the defined terms of
% the A vector. Subsequently the values of the new parameters are defined.
%
% [output] = reducereac(output,indexkeep)
% [output] = reducereac(output,indexkeep,keeporigparameters)

% Information:
% ============
% Systems Biology Toolbox for MATLAB
% Copyright (C) 2005-2006 Henning Schmidt, FCC, henning@sbtoolbox.org
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% check if symbolic toolbox is present
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~symbolicpresent,
    error('The model reduction feature requires the presence of the symbolic toolbox.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Variable input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
keeporigparameters = 0;
if nargin == 3,
    keeporigparameters = varargin{1};
end

output.keeporigparameters = keeporigparameters;
output.reaction_red = [];
output.reaction_red.formula = [];
output.reaction_red.formulaop = [];
output.reaction_red.parameters = [];
output.reaction_red.parametervalues = [];
output.reaction_red.variables = [];
output.reaction_red.variablevalues = [];
output.reaction_red.reactionvalues = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Build new reduced reaction expression ... use new or original parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% First check if new or original parameters (output.keeporigparameters = 0
% does not keep them =1 keeps them eventually, =2 keeps them always).
% copy flag to extra variable to keep it on the initial value in the output
% structure. =1: Original parameters are kept if it means that this leads to
% fewer parameters. 
keeporigparameters = output.keeporigparameters;
if keeporigparameters == 1,
    testparam = {};
    for k = 1:length(indexkeep),
        % determine parameters
        testparam = getNewParameters(testparam,output.reaction_trans.c{indexkeep(k)},output.reaction_orig.parameters);
    end
    if length(testparam) > length(indexkeep),
        % Keeping original parameters would lead to more parameters => skip it
        keeporigparameters = 0;
    end
end

expression = sym(0);
newparam = {};
if keeporigparameters == 0,
    % Use new parameters (eventually more degrees of freedom and a better fit
    newbaseparamname = ['K_' output.reaction '_'];
    for k = 1:length(indexkeep),
        newparam{end+1} = sprintf('%s%d',newbaseparamname,k);
        expression = expression + sym(output.reaction_trans.A{indexkeep(k)})*sym(newparam{end});
    end
else
    % Use original parameters (eventually fewer degrees of freedom and a worse fit)
    for k = 1:length(indexkeep),
        % determine parameters
        newparam = getNewParameters(newparam,output.reaction_trans.c{indexkeep(k)},output.reaction_orig.parameters);
        % build new expression
        expression = expression + sym(output.reaction_trans.A{indexkeep(k)})*sym(output.reaction_trans.c{indexkeep(k)});
    end
end
expression = [output.reaction_trans.b{1} '=' char(expression)];
reacredRHS = solve(expression,output.reaction);
% add info to structure
output.reaction_red.formula = char(reacredRHS);
output.reaction_red.parameters = newparam;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Check if all numerators have been deleted (leads to an error)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if output.reaction_red.formula == '0',
    error('All numerator terms have been deleted from the reaction expression.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the parameter values (starting guesses for eventual optimization)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
if keeporigparameters == 0,
    % New parameters: Use nonlinear regression to estimate the new parameter values
    %Ared = output.reaction_trans.Anum(:,indexkeep);
    %bred = output.reaction_trans.bnum;
    %output.reaction_red.parametervalues = (pinv(Ared)*bred)';
    % Get new parameters from the old ones (as starting guess)
    output.reaction_red.parametervalues = output.reaction_trans.cnum(indexkeep);
else
    % Old parameters: just use the original values
    output.reaction_red.parametervalues = SBparameters(output.model,output.reaction_red.parameters)';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the variables that are present and also get their numeric values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
output.reaction_red.variables = setdiff(explodePCSB(findsym(sym(output.reaction_red.formula))),output.reaction_red.parameters);
variablevalues = [];
for k = 1:length(output.reaction_red.variables),
    index = strmatch(output.reaction_red.variables{k}, output.reaction_orig.variables, 'exact');
    variablevalues(:,end+1) = output.reaction_orig.variablevalues(:,index);
end
output.reaction_red.variablevalues = variablevalues;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the reaction rate for defined experiments 
% for reduced expression
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
% Define all the variables in the workspace
for k = 1:length(output.reaction_red.variables),
    eval(sprintf('%s = output.reaction_red.variablevalues(:,k);',output.reaction_red.variables{k}));
end
% Define the new parameters in the workspace
for k = 1:length(output.reaction_red.parameters),
    eval(sprintf('%s = output.reaction_red.parametervalues(k);',output.reaction_red.parameters{k}));
end
% Evaluate the reduced reaction formula
reaceval = formula2vec(output.reaction_red.formula);
redreactionvalues = eval(reaceval);
output.reaction_red.reactionvalues = redreactionvalues;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine new parameters in expression for parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [newparam] = getNewParameters(newparam, expression, possibleparam)
addparam = explodePCSB(findsym(sym(expression)));
newparam = unique({newparam{:}, addparam{:}});
return
