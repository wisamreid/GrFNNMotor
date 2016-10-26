function varargout = layout(varargin)
% LAYOUT Application M-file for layout.fig
%    FIG = LAYOUT launch layout GUI.
%    LAYOUT('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 09-Jan-2002 12:20:38
global gds driver_window;
if nargin == 0  % LAUNCH GUI
    set(0,'ShowHiddenHandles','on');
    h=findobj('type','figure','Tag','layout_fig');
    delete(h);
    set(0,'ShowHiddenHandles','off');  

	fig = openfig(mfilename,'new');

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
    feval(gds.gui_load_layout,handles);
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
function layoutbox_Callback(h, eventdata, handles, varargin)
global gds ;
index_selected=get(handles.layoutbox,'Value');
list=get(handles.layoutbox,'String');
feval(gds.gui_layoutbox,list,index_selected);
feval(gds.gui_load_layout,handles);
  

% --------------------------------------------------------------------
function layoutOK_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.layoutOK.
delete(handles.layout_fig);
numeric;
