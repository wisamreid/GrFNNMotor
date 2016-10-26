function varargout = pausecont(varargin)
% PAUSECONT Application M-file for pausecont.fig
%    FIG = PAUSECONT launch pausecont GUI.
%    PAUSECONT('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 30-Jan-2002 16:47:58
global sys driver_window;

if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
    loadpausecont(handles);

	% Wait for callbacks to run and window to be dismissed:
	uiwait(fig);

	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
		errordlg(lasterr,mfilename);
        delete(driver_window);
	end

end





% --------------------------------------------------------------------
function special_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.special.
global sys;
val = get(handles.special,'Value');
if val==1
    set(handles.each_point,'Value',0);
    sys.gui.pauseeachpoint = 0;
    set(handles.never,'Value',0);
    sys.gui.pausenever     = 0;
end
sys.gui.pausespecial = val;

% --------------------------------------------------------------------
function each_point_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.each_point.
global sys;
val = get(handles.each_point,'Value');
if val==1
    set(handles.special,'Value',0);
    sys.gui.pausespecial = 0;
    set(handles.never,'Value',0);
    sys.gui.pausenever   = 0;
end
sys.gui.pauseeachpoint = val;

% --------------------------------------------------------------------
function never_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.never.
global sys;
val = get(handles.never,'Value');
if val==1
    set(handles.each_point,'Value',0);
    sys.gui.pauseeachpoint = 0;
    set(handles.special,'Value',0);
    sys.gui.pausespecial   = 0;
end
sys.gui.pausenever = val;

% --------------------------------------------------------------------
function OKpause_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.OKpause.
delete(handles.pausecont)

%---------------------------------------------------------------------
function loadpausecont(handles)
global sys
set(handles.each_point,'Value',sys.gui.pauseeachpoint);
set(handles.never,'Value',sys.gui.pausenever);
set(handles.special,'Value',sys.gui.pausespecial);