function [exp] = SBexp(varargin)
% SBexp: creates an experiment object defining experiment settings 
%
% USAGE:
% ======
% [exp] = SBexp()               creates an empty SBexp object
% [exp] = SBexp(SBstructure)    creates an SBexp object from a MATLAB
%                               structure in the internal experiment format
% [exp] = SBexp(expin)          construction from a given SBexp object (expin)
% [exp] = SBexp('file.exp')     converting a experiment text description 
%                               to an SBexp object.
%
% Output Arguments:
% =================
% exp: SBexp object 

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
    if strcmp('SBexp',class(varargin{1})),
        inputType = 'SBexp';
        expInput = varargin{1};
    elseif strcmp('struct',class(varargin{1})),
        inputType = 'SBstructure';
        SBstructure = varargin{1};
    elseif strcmp('char',class(varargin{1})),
        % check if '.txt' given as extension. If yes, then import text
        % description
        filename = varargin{1};
        if ~isempty(strfind(filename,'.exp')),
            inputType = 'TextExpFile';
        else
            error('Input argument of unknown type');
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
    % parameter settings substructure
    parametersettingsStruct = struct('name',{},'formula',{},'notes',{});
    % initial conditions substructure
    initialconditionsStruct = struct('name',{},'formula',{},'notes',{});    
    % parameter settings substructure
    parameterchangesStruct = struct('name',{},'formula',{},'notes',{});
    % event assignment substructure
    eventassignmentStruct = struct('variable',{},'formula',{});
    % state events substructure
    stateeventsStruct = struct('name',{},'trigger',{},'assignment',eventassignmentStruct,'notes',{});
    % Create SBstructure
    SBstructure = struct('name','unnamed_experiment','notes','no notes','parametersettings',parametersettingsStruct,'initialconditions',initialconditionsStruct,'parameterchanges',parameterchangesStruct,'stateevents',stateeventsStruct);
    % construct the model object
    exp = class(SBstructure,'SBexp');
elseif strcmp('SBexp',inputType),
    % copy the model object
    exp = expInput;
elseif strcmp('SBstructure',inputType),
    % check if the given structure is a SBstructure (only check the
    % top-level fields)
    checkFields = {'name','notes','parametersettings','initialconditions','parameterchanges','stateevents'};
    for k = 1:length(checkFields),
        if ~isfield(SBstructure,checkFields{k}),
            error('Given structure is not a valid internal SBexp structure.');
        end
    end
    % construct the model object
    exp = class(SBstructure,'SBexp');
elseif strcmp('TextExpFile',inputType),
    % check if a file with given filename exists
    if ~exist(filename),
        errorMsg = sprintf('Text model file does not exist in current directory.');
        error(errorMsg);
    end
    % If file exists then first load it
    fid = fopen(strcat(filename),'r');
    expText = fread(fid,inf,'uint8=>char')';
    fclose(fid);
    % then convert it to SBstructure
    [SBstructure, errorMsg] = convertTextToExpSB(expText);
    % Check if error occurred while importing the SBML model
    if length(errorMsg) ~= 0,
        error(errorMsg);
    end
    % construct the model object
    exp = class(SBstructure,'SBexp');
else
    error('Wrong input arguments.');
end
return
