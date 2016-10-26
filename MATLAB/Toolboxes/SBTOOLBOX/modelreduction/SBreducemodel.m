function [reducedmodel] = SBreducemodel(varargin)
% SBreducemodel: Reduces an SBmodel by identifying algebraic relations 
% between time dependent variable, defined by differential equations, and 
% deleting the dependent variables. This function first checks if the 
% SBmodel contains algebraic relations between the state variables. If 
% there are, it lets the user decide which variables to choose as dependent
% ones and to replace the ODEs by algebraic relations.
%
% Due to moiety conservations algebraic relations are often present in 
% biological systems. However, for several analysis methods models are 
% required to be non-singular. This function helps to avoid this problem.
%
% USAGE:
% ======
% [reducedmodel] = SBreducemodel(model)         
% [reducedmodel] = SBreducemodel(model,tol)         
%
% model: SBmodel model to reduce
% tol: tolerance to be used to determine moiety conservations
%
% DEFAULT VALUES:
% ===============
% tol: same as for SBmoietyconservations
%
% Output Arguments:
% =================
% reducedmodel: Reduced SBmodel. If no algebraic relations were present in
%       the model, the initial model is returned unchanged.

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% NEEDED GLOBAL VARIABLES
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global ODEfctname 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK IF SBMODEL OR FILENAME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp('SBmodel',class(varargin{1})),
    % SBmodel
    sbm = varargin{1};
    % Create temporary ODE file
    [ODEfctname, ODEfilefullpath] = SBcreateTempODEfile(sbm);    
else
    error('This function can not be applied to ODE files.');
end

if nargin == 1,
    tol = 0;
elseif nargin == 2,
    tol = varargin{2};
else
    error('Incorrect number of input arguments.');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK AND REDUCE THE MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
reducedmodel = reduceModel(sbm,tol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DELETE FILE 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deleteTempODEfileSB(ODEfilefullpath);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK AND REDUCE THE MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [sbmReduced] = reduceModel(sbm,tol)
% needed global variables
global ODEfctname

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DETERMINE THE MOIETY CONSERVATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
initialCondition = feval(ODEfctname);
[depVarIndex, depVarConstant, depVarFactor, message] = SBmoietyconservations(sbm,SBinitialconditions(sbm),tol);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK IF REDUCTION NECESSARY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(depVarIndex),
    % System is non-singular - return the original system
    disp('Model is non-singular. Nothing is done.');
    sbmReduced = sbm; 
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK THE CONSERVATION EQUATIONS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To be able to include them as variables in an SBmodel it is important
% that the moietyconservations are independent of each other. No algebraic
% loop is allowed.
% for k1 = 1:size(depVarFactor,1),
%     for k2 = 1:size(depVarFactor,1),
%         if k1 ~= k2,
%            if ~isempty(intersect(find(depVarFactor(k1,:)~=0),find(depVarFactor(k2,:)~=0))),
%                error(sprintf('The model can not be reduced due to dependency among the conservations.\nIn the future the toolbox will support these kind of things, but not today.\nSorry!'));
%            end
%         end
%     end
% end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONVERSION OF THE MOIETY CONSERVATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% need to convert these results to a different format.
% result format: dependentVariable = constantTerm + factor1*variable1 + ... + factorn*variablen
% needed format: constantTerm = dependentVariable + NewFactor1*variable1 -
% ... + NewFactorn*variablen
% => change sign in depVarFactor and add a 1 in the places corresponding to
% the dependent variables.
depVarFactor = -depVarFactor;
for k = 1:length(depVarIndex),
    depVarFactor(k,depVarIndex(k)) = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PROCESS SYSTEM
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% System is singular - process the algebraic relations
d = length(depVarIndex);
if d == 1,
    disp(sprintf('The model is singular. %d state needs to be removed to obtain a non-singular system.',d));
else
    disp(sprintf('The model is singular. %d states need to be removed to obtain a non-singular system.',d));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY ALL ALGEBRAIC RELATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(message);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DECIDE WHICH STATES TO REMOVE (USER)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
text = sprintf('You need to decide, for each algebraic relation, which state to remove from the model.\n');
disp(text);
variableNames = feval(ODEfctname,'states');
statesRemove = [];
% cycle through all the algebraic relations, display them and let the
% user decide which state to remove.
for k1 = 1:length(depVarConstant),
    indices = setdiff(find(depVarFactor(k1,:) ~= 0),statesRemove);
    text = sprintf('Relation %d: %s',k1,getTextAlgebraicRelation(depVarConstant(k1),depVarFactor(k1,:)));
    disp(text);
    text = sprintf('Choose the number of the state that should be removed (1-%d).',length(indices));
    disp(text);
    text = '';
    for k2 = 1:length(indices),
        text = sprintf('%s (%d)%s  ',text,k2,variableNames{indices(k2)});
    end
    disp(text);
    indexRemove = input('State to remove: ');
    if indexRemove > length(indices) || indexRemove < 1,
        error('Wrong input!');
    end
    statesRemove = [statesRemove indices(indexRemove)];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCT REDUCED MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% get the structure of the SBmodel
modelStructure = SBstruct(sbm); % will contain the reduced information
fullModelStructure = SBstruct(sbm); % contains the full information (needed only for updating the additional information fields of the variables)
% create new "states" entry and copy only the states that are not
% deleted into it
states = [];
statesCount = 1;
for k1 = 1:length(modelStructure.states),
    if sum(ismember(statesRemove,k1)) == 0,
        % copy this state information to new model
        states(statesCount).name = modelStructure.states(k1).name;
        states(statesCount).initialCondition = modelStructure.states(k1).initialCondition;
        states(statesCount).ODE = modelStructure.states(k1).ODE;
        states(statesCount).type = modelStructure.states(k1).type;
        states(statesCount).compartment = modelStructure.states(k1).compartment;
        states(statesCount).unittype = modelStructure.states(k1).unittype;
        states(statesCount).notes = modelStructure.states(k1).notes;
        statesCount = statesCount + 1;
    end
end
% Include the state information in the model structure
modelStructure.states = states;
% Now add the algebraic relations as variables.
% It is important to add the relations before all other variables!
% Since they might be used in the other variables!
variables = [];
variableCount = 1;
parameterCount = length(modelStructure.parameters)+1;
% First add the algebraic relations
for k1 = 1:length(depVarConstant),
    % Construct variable formula corresponding to algebraic relation
    % Get the constant term
    constantTerm = depVarConstant(k1);
    % Scale the factor vector with -1 and set reduced element to zero
    factor = depVarFactor(k1,:);
    scaleFactor = 1/abs(factor(statesRemove(k1)));
    factor = -factor*scaleFactor;
    factor(statesRemove(k1)) = 0;
    % Scale the constant term
    constantTerm = constantTerm*scaleFactor;
    % get the variableName
    variable = variableNames{statesRemove(k1)};
    % Construct the formula string
    moietyConservationName = sprintf('conserv%d',k1);
    formula = sprintf('%s',moietyConservationName);
    for k2 = 1:length(factor),
        if factor(k2) > 0,
            % here 0 can be used (zerocheck already done)
            formula = sprintf('%s+%g*%s',formula,abs(factor(k2)),variableNames{k2});
        elseif factor(k2) < 0,
            % here 0 can be used (zerocheck already done)
            formula = sprintf('%s-%g*%s',formula,abs(factor(k2)),variableNames{k2});
        end
    end
    % Add algebraic relation as variable
    variables(variableCount).name = variable;
    variables(variableCount).formula = formula;
    variables(variableCount).notes = 'Removed algebraic relation';
    variables(variableCount).type = fullModelStructure.states(statesRemove(k1)).type;
    variables(variableCount).compartment = fullModelStructure.states(statesRemove(k1)).compartment;
    variables(variableCount).unittype = fullModelStructure.states(statesRemove(k1)).unittype;
    variableCount = variableCount + 1;
    % add the constant term as a parameter
    modelStructure.parameters(parameterCount).name = moietyConservationName;
    modelStructure.parameters(parameterCount).value = str2num(sprintf('%g',constantTerm));
    moietyConservationNotes = sprintf('Moiety conservation for: %s',getTextAlgebraicRelation(depVarConstant(variableCount-1),depVarFactor(variableCount-1,:)));
    moietyConservationNotes = moietyConservationNotes(1:strfind(moietyConservationNotes,'=')-1);
    modelStructure.parameters(parameterCount).type = 'isParameter';
    modelStructure.parameters(parameterCount).compartment = '';
    modelStructure.parameters(parameterCount).unittype = '';
    modelStructure.parameters(parameterCount).notes = moietyConservationNotes;
    parameterCount = parameterCount + 1;
end
% Now add the other variables
for k = 1:length(fullModelStructure.variables),
    variables(variableCount).name = fullModelStructure.variables(k).name;
    variables(variableCount).formula = fullModelStructure.variables(k).formula;
    variables(variableCount).type = fullModelStructure.variables(k).type;
    variables(variableCount).compartment = fullModelStructure.variables(k).compartment;
    variables(variableCount).unittype = fullModelStructure.variables(k).unittype;
    variables(variableCount).notes = fullModelStructure.variables(k).notes;
    variableCount = variableCount + 1;
end
% Include the variable information in the model structure
modelStructure.variables = variables;
% Create SBmodel with reduced model structure
sbmReduced = SBmodel(modelStructure);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DISPLAY THE ALGEBRAIC RELATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [text] = getTextAlgebraicRelation(depVarConstant,depVarFactor)
    % needed global variables
    global ODEfctname
    % get variable names from ODE file
    variableNames = feval(ODEfctname,'states');
    text = sprintf(' = %g',depVarConstant);
    for k2 = 1:length(depVarFactor),
        if abs(depVarFactor(k2)) > 0,    % here 0 can be used (zerocheck already done)
            if depVarFactor(k2) > 0
                element = sprintf('+ %g',depVarFactor(k2));
            else 
                element = sprintf('- %g',abs(depVarFactor(k2)));
            end
            text = sprintf('%s %s %s',element,variableNames{k2},text);
        end
    end
return

