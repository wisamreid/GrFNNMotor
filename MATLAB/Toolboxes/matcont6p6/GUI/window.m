function varargout = window(varargin)
% WINDOW Application M-file for window.fig
%    FIG = WINDOW launch window GUI.
%    WINDOW('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 12-Dec-2002 14:45:27
global np sys MC;
if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
    guidata(fig, handles);
    set(handles.npoints,'String',num2str(sys.gui.plot_points));
    np=sys.gui.plot_points;
	% Wait for callbacks to run and window to be dismissed:
	uiwait(fig);

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
function npoints_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.npoints
npoints=str2double(get(handles.npoints,'String'));
if npoints<=0
    warndlg('The number must be positive and non-zero');
    set(handles.npoints,'String',1);
    return
else
    set(handles.npoints,'String',num2str(round(npoints)));
end

% --------------------------------------------------------------------
function OK_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.OK.
global sys path_sys MC;
%set(0,'ShowHiddenHandles','on');
%f=findobj('Type','Uimenu','Tag','options_window');
%set(0,'ShowHiddenHandles','off');
sys.gui.plot_points=str2double(get(handles.npoints,'String'));
set(MC.mainwindow.options_window,'Userdata',sys.gui.plot_points);
delete(handles.window);

% --------------------------------------------------------------------
function Cancel_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.Cancel.
global np sys path_sys MC;
sys.gui.plot_points=np;
%set(0,'ShowHiddenHandles','on');
%f=findobj('Type','uimenu','tag','options_window');
%set(0,'ShowHiddenHandles','off');
set(MC.mainwindow.options_window,'Userdata',np);
delete(handles.window);
