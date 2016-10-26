function varargout = archive(varargin)
% ARCHIVE Application M-file for archive.fig
%    FIG = ARCHIVE launch archive GUI.
%    ARCHIVE('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 13-Dec-2002 15:59:09
global sys np
if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
    set(handles.edit1,'String',num2str(sys.gui.ncurves));
    np=sys.gui.ncurves;
	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
		disp(lasterr);
	end

end





% --------------------------------------------------------------------
function OK_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.OK.
global sys
sys.gui.ncurves=str2double(get(handles.edit1,'String'));
delete(handles.archive);

% --------------------------------------------------------------------
function Cancel_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.Cancel.
global sys np
sys.gui.ncurves=np;
delete(handles.archive);

% --------------------------------------------------------------------
function edit1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit1.
ncurves=str2double(get(handles.edit1,'String'));
if ncurves<=0
    warndlg('The number must be positive and non-zero');
    set(handles.edit1,'String',2);
    return
else
    set(handles.edit1,'String',num2str(round(ncurves)));
end
