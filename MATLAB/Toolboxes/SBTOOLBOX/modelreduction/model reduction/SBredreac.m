function [modelred] = redreac(output,reaction,varargin)
% redreac: Reduction of single reaction expressions with the goal of reducing
% complexity and number of parameters. This function requires the prior
% execution of the prepredreac.
% 
% USAGE:
% ======
% [modelred] = redreac(output,reaction)
% [modelred] = redreac(output,reaction,OPTIONS)
%
% output: Result returned from prepredreac
% reaction: Name of the reaction that should be reduced
% OPTIONS: structure containing options for the algorithm:
%        OPTIONS.tol: tolerance for singularity detection (smalles SV)
%        OPTIONS.keeporigparameters: =0: do not keep original parameters
%                                    =2: do always keep original parameters
%                                    =1: keep original parameters only if
%                                        it leads to fewer parameters
%        OPTIONS.numeratorweighting: =1: weight numerator terms and denumerator terms 
%                                    such that numerators are kept
%                                    =0: dont do any weighting
%
% DEFAULT VALUES:
% ===============
% OPTIONS.tol:                  1e-6
% OPTIONS.keeporigparameters:   0
% OPTIONS.numeratorweighting:   0
%
% Output Arguments:
% =================
% modelred: SBmodel in which the defined reaction is reduced

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Handle variable input arguments (options)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OPTIONS = [];
if nargin == 3,
    OPTIONS = varargin{1};
end
% default values
tol = 1e-6;              % tolerance for singularity detection (smalles SV)
keeporigparameters = 0;  % do not keep original parameters
numeratorweighting = 0;  % do not weight numerators/denumerators depending on their number
% tol
if isfield(OPTIONS,'tol'),
    if ~isempty(OPTIONS.tol),
        tol = OPTIONS.tol;
    end
end
% keeporigparameters
if isfield(OPTIONS,'keeporigparameters'),
    if ~isempty(OPTIONS.keeporigparameters),
        keeporigparameters = OPTIONS.keeporigparameters;
    end
end
% numeratorweighting
if isfield(OPTIONS,'numeratorweighting'),
    if ~isempty(OPTIONS.numeratorweighting),
        numeratorweighting = OPTIONS.numeratorweighting;
    end
end

disp(' ');
disp('#####################################################################');
disp(sprintf('# Reduction of reaction:  %s', reaction));
disp('#####################################################################');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine reduction information for specified reaction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
output2 = reducereactionprep(output,reaction);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine which terms to reduce
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% number terms in original equation
nrorig = size(output2.reaction_trans.A,2);
indexkeep = [1:nrorig];
indexkeep_backup = {};
while 1,
    try
        % Reduce (do this also for the case where indexkeep contains all terms)
        output3 = reducereac(output2,indexkeep,keeporigparameters);
    catch
        % catch errors that can occurr in cases when no reduction is possible
        % example: reaction only defined by constant rate
        %          ...
        
        disp('************************** All nominators might have been cancelled ... try nominator weighting!');
    end
    % Optimize parameter values
    output4 = reduceoptim(output3);
    % check if optimization successful (do it by looking at it!)
    orig = output3.reaction_orig.reactionvalues;
    redopt = output4.reaction_opt.reactionvalues;
    figH = figure(1); clf; 
    plot(orig); hold on; plot(redopt,'r--'); 
    % plot some info
disp(' ');
if length(indexkeep) > 1,
    disp('Identifiabilty and sensitivity measures ');
    disp('#######################################');
    [Us,Ss,Vs] = svd(output3.reductioninfo.Anumscaled(:,indexkeep));
    disp(sprintf('Singular values Mscaled: %s', sprintf('%g ',diag(Ss))));
    disp(sprintf('Elements v*(Mscaled): %s',sprintf('%g ',Vs(:,end))));
    disp(' ');
    [Uc,Sc,Vc] = svd(output3.reductioninfo.AnumC(:,indexkeep));
    disp(sprintf('Singular values Mc: %s', sprintf('%g ',diag(Sc))));
    disp(sprintf('Elements v*(Mc): %s',sprintf('%g ',Vc(:,end))));
    disp(' ');
end
disp('Statistics and preliminary result ');
disp('#################################');
    text = sprintf('Nr parameters in original/reduced expression: %d/%d\n',length(output2.reaction_orig.parameters),length(output3.reaction_red.parameters));
    indexdelete = setdiff([1:nrorig],indexkeep);
    deleteText = '';
    for k = 1:length(indexdelete),
        deleteText = strcat(deleteText, ', ', output3.reaction_trans.A{indexdelete(k)});
    end
    text = sprintf('%sDeleted terms: %s\n',text,deleteText(2:end));
    text = sprintf('%s%s_red = %s\n',text,output2.reaction,output3.reaction_red.formula);
    disp(text);
    % Decide what to do next
    if length(indexkeep_backup) == 0,
        while 1,
            check = input('Continue (2), or end (0) (switch num weighting: 99): ');
            if check == 2 || check == 0,
                break;
            elseif check == 99,
                numeratorweighting = ~numeratorweighting;
                disp(sprintf('numeratorweighting = %d',numeratorweighting));
            end
        end
    else
        while 1,
            check = input('Continue (2), go back (1), end (0) (switch num weighting: 99): ');
            if check == 2 || check == 1 || check == 0,
                break;
            elseif check == 99,
                numeratorweighting = ~numeratorweighting;
                disp(sprintf('numeratorweighting = %d',numeratorweighting));
            end
        end
    end
    if check == 1,
        % go a step back
        indexkeep = indexkeep_backup{end};
        indexkeep_backup = indexkeep_backup(1:end-1);
    elseif check == 0,
        % end reduction
        break;
    else
        % continue
        indexkeep_backup{end+1} = indexkeep;
        % Select what to reduce
        Anumscaled = output2.reductioninfo.Anumscaled(:,indexkeep);
        AnumC = output2.reductioninfo.AnumC(:,indexkeep);
        [U,S,V] = svd(Anumscaled);
        SAns = diag(S)';
        VAns = V(:,end);
        [U,S,V] = svd(AnumC);
        SAnc = diag(S)';
        VAnc = V(:,end);
        V1 = abs(VAns);
        V2 = abs(VAnc);
        % do numeratorweighting (if on then set numerator elements to 0),
        % avoidng the deletion of numerators
        if numeratorweighting == 1,
            isnomterm = output2.reductioninfo.isnomterm(indexkeep);
            numbernom = sum(isnomterm);
            numberden = length(isnomterm)-numbernom;
            factor = [ones(numberden,1); 0*ones(numbernom,1)];
            V1 = V1.*factor;
            V2 = V2.*factor;
        end
        V1 = max(V1-mean(V1),0);
        V2 = max(V2-mean(V2),0);
        if SAns(end) < tol,
            disp('########################################');
            disp('# Reduction:  Based on identifiability #');
            disp('########################################');
            V = V1+0.01*V2;
        else
            disp('##################################');
            disp('# Reduction based on sensitivity #');
            disp('##################################');
            V = V2;
        end
        [dummy,indexdeleteNEXT] = max(V);
        indexkeep(indexdeleteNEXT) = [];   % delete element
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Update the model with the reduced reaction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
modelred= redupdate(output4);
