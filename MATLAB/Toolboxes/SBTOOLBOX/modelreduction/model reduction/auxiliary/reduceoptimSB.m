function [output] = reduceoptim(output)
% reduceoptim: optimizes reduced parameter values

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
% Optimize new parameters with estimates as starting guess
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
global opt_ref opt_formula opt_param opt_var opt_varvalues reac
opt_ref = output.reaction_orig.reactionvalues;           % reaction value from unreduced reaction
opt_formula = output.reaction_red.formula;               % reduced formula
opt_param = output.reaction_red.parameters;              % parameternames in reduced formula
opt_paramvalues = output.reaction_red.parametervalues;   % parameternames in reduced formula
opt_var = output.reaction_red.variables;                 % variable names
opt_varvalues = output.reaction_red.variablevalues;      % variable values
options.maxiter = 100;
options.lowbounds = -1e10*ones(1,length(opt_param));
options.highbounds = 1e10*ones(1,length(opt_param));
options.silent = 1;
old = inf;
FVAL = -inf;
while abs(old-FVAL) > 1e-4
    old = FVAL;
    [opt_paramvalues, FVAL] = simplexSB(@costfunction,opt_paramvalues,options);
end
output.reaction_opt.parametervalues = opt_paramvalues;
output.reaction_opt.reactionvalues = reac;
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
function [cost] = costfunction(X)
global opt_ref opt_formula opt_param opt_var opt_varvalues reac
% assign parameter values
for k = 1:length(X),
    eval(sprintf('%s = %f;',opt_param{k},X(k)));
end
% assign variable values
for k = 1:length(opt_var),
    eval(sprintf('%s = opt_varvalues(:,k);',opt_var{k}));
end
% evaluate reaction equation
reac = eval(formula2vec(opt_formula));
% determine cost
cost = sqrt(sum((reac-opt_ref).^2));
return