function [model] = SBmodel(varargin)
% SBmodel: creates a model object for a biological system
%
% USAGE:
% ======
% [model] = SBmodel()               creates an empty SBmodel 
% [model] = SBmodel(SBstructure)    creates an SBmodel from a MATLAB
%                                   structure in the internal model format
% [model] = SBmodel(modelin)        construction from a given SBmodel (modelin)
% [model] = SBmodel('file.xml')     converting an SBML model file to an 
%                                   SBmodel.
% [model] = SBmodel('file.xml',flag) converting an SBML model file to an 
%                                    SBmodel. 
% flag = 0: standard behavior.
% flag = 1: special case of SBmodel. This functionality is mainly thought
% for the import of (in)completely defined CellDesigner SBML models. It
% does not require that a model is fully defined in terms of parameter
% values, and rate equations, etc. Furthermore it does a name to id
% conversion. This is necessary, since in CellDesigner the ids of the
% elements are chosen automatically and in the SBT the ids are used as
% names for states, variables, etc. It is made sure that all new ids are
% unique. Only for Level 2 SBML models. It might also be used for models
% defined using other modelling environments where the ids have no
% interpretation.
% flag = 2: same as flag 1, additionally all extra SBML related information
% in the SBmodel is taken away.
%
% [model] = SBmodel('textfile.txt') converting a ODE text file description of a 
%                                   model to an SBmodel. The textfile
%                                   description has the format used in
%                                   SBedit. If using this kind of model
%                                   description the extension '.txt' needs to
%                                   be specified!
%
% [model] = SBmodel('textfile.txtbc') converting a text file using a 
%                                   biochemical description of a 
%                                   model to an SBmodel. The textfile
%                                   description has the format used in
%                                   SBeditBC. If using this kind of model
%                                   description the extension '.txtbc' needs to
%                                   be specified!
%
% Output Arguments:
% =================
% model: SBmodel object 

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

flag = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CHECK THE TYPE OF THE INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==0,
    inputType = 'empty';
elseif nargin == 1 || nargin == 2,
    if strcmp('SBmodel',class(varargin{1})),
        inputType = 'SBmodel';
        sbmInput = varargin{1};
    elseif strcmp('struct',class(varargin{1})),
        inputType = 'SBstructure';
        SBstructure = varargin{1};
    elseif strcmp('char',class(varargin{1})),
        % check if '.txt' or '.txtbc' given as extension. If yes, then import text
        % model, otherwise assume an SBML model is to be imported.
        filename = varargin{1};
        if ~isempty(strfind(filename,'.txtbc')),
            inputType = 'TextBCModelFile';
        elseif ~isempty(strfind(filename,'.txt')),
            inputType = 'TextModelFile';
        else
            inputType = 'SBMLmodelFile';
            if nargin == 2,
                flag = varargin{2};
            end
        end
    else 
        error('Input argument of unknown type');
    end
else
    error('Wrong number of input arguments.');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CONSTRUCT THE SBMODEL OBJECT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp('empty',inputType),
    % Create empty SBstructure
    % functions substructure
    functionsStruct = struct('name',{},'arguments',{},'formula',{},'notes',{});
    % states substructure
    statesStruct = struct('name',{},'initialCondition',{},'ODE',{},'type',{},'compartment',{},'unittype',{},'notes',{});
    % parameters substructure
    parametersStruct = struct('name',{},'value',{},'type',{},'compartment',{},'unittype',{},'notes',{});
    % variables substructure
    variablesStruct = struct('name',{},'formula',{},'type',{},'compartment',{},'unittype',{},'notes',{});
    % reactions substructure
    reactionsStruct = struct('name',{},'formula',{},'notes',{},'reversible',{});
    % event assignment substructure
    eventassignmentStruct = struct('variable',{},'formula',{});
    % event substructure
    eventStruct = struct('name',{},'trigger',{},'assignment',eventassignmentStruct,'notes',{});
    % Create SBstructure
    SBstructure = struct('name','unnamed_model','notes','','functions',functionsStruct,'states',statesStruct,'parameters',parametersStruct,'variables',variablesStruct,'reactions',reactionsStruct,'events',eventStruct,'functionsMATLAB','');
    % construct the model object
    model = class(SBstructure,'SBmodel');
elseif strcmp('SBmodel',inputType),
    % copy the model object
    model = sbmInput;
elseif strcmp('SBstructure',inputType),
    % check if the given structure is a SBstructure (only check the
    % top-level fields)
    checkFields = {'name','notes','functions','states','parameters','variables','reactions','events','functionsMATLAB'};
    for k = 1:length(checkFields),
        if ~isfield(SBstructure,checkFields{k}),
            errorMsg = sprintf('Given structure is not a valid internal SBmodel structure.');
            error(errorMsg);
        end
    end
    % construct the model object
    model = class(SBstructure,'SBmodel');
elseif strcmp('SBMLmodelFile',inputType),
    % check if a file with given filename exists
    [path,filename,ext,versn] = fileparts(filename);
    filename = fullfile(path, [filename '.xml']); 
    if ~exist(filename),
        errorMsg = sprintf('SBML model, "%s", does not exist.', filename);
        error(errorMsg);
    end
    % If file exists then import it
    if flag == 0,
        [SBstructure,errorMsg] = SBimportSBML(filename);
        if length(errorMsg) ~= 0, error(errorMsg);  end
    elseif flag == 1 || flag == 2,
        [SBstructure,errorMsg] = SBimportSBMLCD(filename);
        if length(errorMsg) ~= 0, warning(errorMsg);  end
    end
    if flag == 2,
        for k = 1:length(SBstructure.states),
            SBstructure.states(k).type = '';
            SBstructure.states(k).compartment = '';
            SBstructure.states(k).unittype = '';
        end
        for k = 1:length(SBstructure.variables),
            SBstructure.variables(k).type = '';
            SBstructure.variables(k).compartment = '';
            SBstructure.variables(k).unittype = '';
        end
        for k = 1:length(SBstructure.parameters),
            SBstructure.parameters(k).type = '';
            SBstructure.parameters(k).compartment = '';
            SBstructure.parameters(k).unittype = '';
        end
        for k = 1:length(SBstructure.reactions),
            SBstructure.reactions(k).reversible = 0;
        end
    end
    % construct the model object
    model = class(SBstructure,'SBmodel');
elseif strcmp('TextModelFile',inputType),
    % check if a file with given filename exists
    if ~exist(filename),
        errorMsg = sprintf('Text model file does not exist in current directory.');
        error(errorMsg);
    end
    % If file exists then first load it
    fid = fopen(strcat(filename),'r');
    modelText = fread(fid,inf,'uint8=>char')';
    fclose(fid);
    % then convert it to SBstructure
    [SBstructure, errorMsg] = convertTextToModelSB(modelText);
    % Check if error occurred while importing the SBML model
    if length(errorMsg) ~= 0,
        error(errorMsg);
    end
    % construct the model object
    model = class(SBstructure,'SBmodel');
elseif strcmp('TextBCModelFile',inputType),
    % check if a file with given filename exists
    if ~exist(filename),
        errorMsg = sprintf('TextBC model file does not exist in current directory.');
        error(errorMsg);
    end
    % If file exists then first load it
    fid = fopen(strcat(filename),'r');
    modelText = fread(fid,inf,'uint8=>char')';
    fclose(fid);
    % then convert it to SBstructure
    [SBstructure, errorMsg] = convertTextBCToModelSB(modelText);
    % Check if error occurred while importing the SBML model
    if length(errorMsg) ~= 0,
        error(errorMsg);
    end
    % construct the model object
    model = class(SBstructure,'SBmodel');
else
    error('Wrong input arguments.');
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DELETE WHITESPACES IN STRINGS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Useful for taking away whitespaces in kineticLaw formulas, as
% seen in some example models
function [outputString] = removeWhiteSpace(inputString)
outputString = strrep(inputString,' ','');
% return
return