function [N,stateNames] = reverseMoietyConservationsSB(sbm,varargin)
% reverseMoietyConservationsSB: Function searching for moietyconservations
% in the system. Then an augmented stoichiometric matrix is calculated in 
% which also the species determined by moiety conservations appear.
%
% input arguments: sbm, specieFlag (optional)
% if the flag is set only variables that are species will be considered

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

specieFlag = 0;
if nargin == 2,
    specieFlag = varargin{1};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GET INFORMATION ABOUT EVENTUAL MOIETY CONSERVATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This is done by searching all variables for a variable defined in the
% format of:    VARIABLE = PARAMETER - factor*STATE ... - factor*STATE ... - factor*STATE 		
% Then it is assumed that VARIABLE + factor*STATE + ... is a conserved moiety and
% stored as such for later use. The "factor*" is optional but if existent
% has to be numerical.

if ~specieFlag,
    % just consider all variables
    [allVariables, allVariableFormulas] = SBvariables(sbm);
else
    allVariables = {};
    allVariableFormulas = {};
    sbmstruct = struct(sbm);
    for k=1:length(sbmstruct.variables),
        if strcmp(sbmstruct.variables(k).type,'isSpecie'),
            allVariables{end+1} = sbmstruct.variables(k).name;
            allVariableFormulas{end+1} = sbmstruct.variables(k).formula;
        end
    end
end

moietyconservations = [];
indexMC = 1;
for k1 = 1:length(allVariables),
    variableFormula = allVariableFormulas{k1};
    terms = explodePCSB(variableFormula,'-');
    if length(terms) > 1,
        % process the terms only in the case that there are more than one single term.
        % first check if the first term is composed only of a parameter
        % name or a numeric value
        firstTerm = strtrim(terms{1});
        if isparameterSB(sbm,firstTerm) || ~isempty(str2num(firstTerm)),
            % only continue if first term was a parameter
            % now cycle through the remaining terms and check if they 
            % are states or factor*states
            formatCorrect = 1;
            listofstates = {};  % list of states in conservation
            listoffactors = []; % list of factors in conservation
            indexstates = 1;
            for k2 = 2:length(terms),
                checkStateTerms = explodePCSB(terms{k2},'*');
                % only accept lengths 1 or 2 (possibly with or without
                % factor) otherwise break
                if length(checkStateTerms) == 1,
                    if isstateSB(sbm,strtrim(checkStateTerms{1})),
                        % yes, state found => add it to the list
                        listofstates{indexstates} = strtrim(checkStateTerms{1});
                        listoffactors(indexstates) = 1;
                        indexstates = indexstates + 1;
                    else
                        formatCorrect = 0;
                        break;
                    end
                elseif length(checkStateTerms) == 2,
                    % second term needs to be a state
                    % if this is the case the first term needs to be
                    % numeric
                    if isstateSB(sbm,strtrim(checkStateTerms{2})),
                        % yes, state found => check if first term is
                        % numeric
                        value = str2num(checkStateTerms{1});
                        if ~isempty(value),
                            listofstates{indexstates} = strtrim(checkStateTerms{2});
                            listoffactors(indexstates) = value;
                            indexstates = indexstates + 1;
                        else
                            formatCorrect = 0;
                            break;
                        end
                    else
                        formatCorrect = 0;
                        break;
                    end
                else
                    formatCorrect = 0;
                    break;
                end
            end
            % if formatCorrect = 1 use it as moiety conservation
            if formatCorrect == 1,
                % it has the right format for defining a moiety
                % conservation => accept it
                moietyconservations(indexMC).variablename = allVariables{k1};
                moietyconservations(indexMC).states = listofstates;
                moietyconservations(indexMC).factors = listoffactors;
                indexMC = indexMC + 1;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONVERT A NONSINGULAR SYSTEM TO A SINGULAR IN CASE OF MOIETY CONSERVATIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% first get the stoichiometric matrix of the system
% Determine the stoichiometric matrix and the state names
[N,stateNames] = SBstoichiometry(sbm,1,1);
% now differentiate the moietyconservations and add the components defined
% by variables to the stoichiometric matrix
for k1 = 1:length(moietyconservations),
    % construct the row for the stoichiometric matrix, belonging to the variable
    Nrownew = zeros(1,length(SBreactions(sbm)));
    for k2 = 1:length(moietyconservations(k1).states),
        name = moietyconservations(k1).states{k2};
        % find the index of the statename
        stateIndex = strmatch(name,stateNames,'exact');
        % get its row of the stoichiometric matrix multiplied by the factor
        % of the state in the moiety conservation and substract from
        % Nrownew
        Nrownew = Nrownew - moietyconservations(k1).factors(k2)*N(stateIndex,:);
    end
    % add species to stateNames (its not a state in the original model but
    % in the singular version that is constructed now)
    stateNames{end+1} = moietyconservations(k1).variablename;
    % add new row to N
    N = [N; Nrownew];
end

return