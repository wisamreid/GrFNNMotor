function [MAstructure] = SBconvert2MA(model)
% SBconvert2MA: Takes an SBmodel and if possible returns the stoichiometric
% matrix, the numeric kinetic parameters, and the initial conditions. The
% SBmodel is not allowed to contain variables, functions, events, or
% functionsMATLAB. All right hand sides of the ODEs have to be defined in
% terms of reaction names. All reactions need to be irreversible and mass action type! 
%
% Please note that the function is not parsing the reaction kinetics and
% checking if they are correct MA kinetic rate laws. Instead, the function
% expects the reaction kinetics given in the following format:
%       kineticParameter*...
% for example:
%       kineticParameter*speciesName1   (speciesName1*kineticParameters is
%                                        not allowed!)      
%       kineticParameter*speciesName1^2
%       kineticParameter*speciesName1*speciesName2
%       etc.
% If no kinetic parameter is present as first element, the function warns the user
% and assumes a default value of 1. 
%
% USAGE:
% ======
% [MAstructure] = SBconvert2MA(model) 
%
% model: SBmodel 
%
% Output Arguments:
% =================
% MAstructure: structure containing information about the MA model
%          MAstructure.N: stoichiometric matrix
%          MAstructure.L: reactant stoichiometric matrix
%          MAstructure.kineticParameters: vector containing the kinetic
%               parameters for each reaction
%          MAstructure.initialConditions: vector of initial conditions for all
%               species
%          MAstructure.species: cell-array with species names
%          MAstructure.reactions: cell-array with reaction names

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Process the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% make all reactions irreversible (if necessary)
modelirr = SBmakeirreversible(model);
% take away parenthesis around backward reaction expressions
ms = struct(modelirr);
for k = 1:length(ms.reactions),
    ms.reactions(k).formula = strrep(ms.reactions(k).formula,'(','');
    ms.reactions(k).formula = strrep(ms.reactions(k).formula,')','');
end
model = SBmodel(ms);
SBstructure = SBstruct(model);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK IF NON DESIRED ELEMENTS ARE PRESENT IN THE MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% undesired: 
%            - functions
%            - variables
%            - events
%            - functionsMATLAB
%            - states defined by rule
% check functions, variables, events, functionsMATLAB
if length(SBstructure.functions) ~= 0,
    error('The model contains functions. This is not allowed when converting it to an MA representation.');
end
if length(SBstructure.variables) ~= 0,
    error('The model contains variables. This is not allowed when converting it to an MA representation.');
end
if length(SBstructure.events) ~= 0,
    error('The model contains events. This is not allowed when converting it to an MA representation.');
end
if length(SBstructure.functionsMATLAB) ~= 0,
    error('The model contains MATLAB functions. This is not allowed when converting it to an MA representation.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK FOR REVERSIBLE REACTIONS (NOT ALLOWED)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[reactionNames, reactionFormulas] = SBreactions(model);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET THE (REACTANT) STOICHIOMETRIC MATRIX 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% also check if states are defined by rules, that is, not by summation of
% reaction terms. easily done by checking the stoichiometric matrix and 
% comparing number of components in it and number of all states in the
% model.
allStates = SBstates(model);
[N,componentNames,reactionNames,reversibleFlag] = SBstoichiometry(model);
if length(allStates) ~= length(componentNames),
    error('Not all ODEs are expressed in terms of reactions. The full stoichiometric matrix can thus not be determined.');
end
L = SBreactantstoichiometry(model);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET THE KINETIC PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[parameterNames] = SBparameters(model);
% evaluate all reaction rates for unity species values. should give the
% same result as the previous method. if not then error!
% use this check only in case that there are as many reactions as
% parameters.
if length(reactionNames) == length(parameterNames),
   [x,y,z,kineticParameters] = SBreactions(model,ones(1,length(allStates)));
else
    [x,y,z,kineticParameters] = SBreactions(model,ones(1,length(allStates)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET THE INITIAL CONDITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initialconditions = SBinitialconditions(model);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FORMAT OUTPUT VALUES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MAstructure.N = N;
MAstructure.L = L;
MAstructure.kineticParameters = kineticParameters(:);
MAstructure.initialConditions = initialconditions(:);
MAstructure.species = allStates;
MAstructure.reactions = {SBstructure.reactions.name}';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RETURN
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return



