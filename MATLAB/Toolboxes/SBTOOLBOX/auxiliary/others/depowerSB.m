function [model,changeFlag] = depowerSB(model)
% depowerSB: simple function that replaces power(x,y) expressions in models
% to (x)^(y)
%
% USAGE:
% ======
% model = depowerSB(model)
%
% Output Arguments:
% =================
% model:        SBmodel with replaced power expressions
% changeFlag:   =0 => no changes, =1 => changes made

global changeFlag
changeFlag = 0;
sbstruct = SBstruct(model);
% states
for k = 1:length(sbstruct.states),
    sbstruct.states(k).ODE = exchangepowerexp(sbstruct.states(k).ODE);
end
% variables
for k = 1:length(sbstruct.variables),
    sbstruct.variables(k).formula = exchangepowerexp(sbstruct.variables(k).formula);
end
% reactions
for k = 1:length(sbstruct.reactions),
    sbstruct.reactions(k).formula = exchangepowerexp(sbstruct.reactions(k).formula);
end
% events
for k = 1:length(sbstruct.events),
    sbstruct.events(k).trigger = exchangepowerexp(sbstruct.events(k).trigger);
    for k2 = 1:length(sbstruct.events(k).assignment),
        sbstruct.events(k).assignment(k2).formula = exchangepowerexp(sbstruct.events(k).assignment(k2).formula);
    end        
end
% functions
for k = 1:length(sbstruct.functions),
    sbstruct.functions(k).formula = exchangepowerexp(sbstruct.functions(k).formula);
end
model = SBmodel(sbstruct);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% regexprep command doing the replacement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [formula] = exchangepowerexp(formula)
global changeFlag
oldformula = formula;
formula = regexprep(['#' formula],'([\W]+)',' $1 ');
formula = regexprep(formula,'[\s]power[\s]*\(([^,]+),([^,]+)\)','($1)^($2)');
formula = regexprep(formula,'\s','');
formula = formula(2:end);
if ~strcmp(oldformula,formula),
    changeFlag = 1;
end
return
    
