function varargout = stop(varargin)
% STOP Application M-file for stop.fig
%    FIG = STOP launch stop GUI.
%    STOP('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.5 12-Jun-2009 09:59:19
global driver_window calculation_progress
if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');
    % Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
    driver_window = fig;
	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
		errordlg(lasterr,mfilename);
        if ishandle(driver_window)
            delete(driver_window);   
        end
	end

end
% --------------------------------------------------------------------
function resume_Callback(h, eventdata, handles, varargin)
global driver_window calculation_progress gds;
    
    set(driver_window,'Userdata',0);
    calculation_progress = 1;

%-------------------------------------------------------------------
function pause_Callback(h, eventdata,handles,varargin)
global driver_window calculation_progress gds;
    set(driver_window,'UserData',1);
    calculation_progress = 2;

%-------------------------------------------------------------------
function stop_Callback(h, eventdata,handles,varargin)
global driver_window calculation_progress gds;
    set(driver_window,'UserData',0);
    calculation_progress = 0;

%-------------------------------------------------------------------
function teststop(h)
global driver_window calculation_progress
    if isempty(driver_window)||calculation_progress==0
        return
    end
    key = double(get(h,'Currentcharacter'));
    if isempty(key),return;end
    switch key
    case 27% press Esc to stop
        stop('stop_Callback',driver_window,[],[]);  
    case 13%press Enter to pause
        stop('pause_Callback',driver_window,[],[]);    
    case 32%press Space bar to resume
        stop('resume_Callback',driver_window,[],[]);    
    end


% --- Executes during object creation, after setting all properties.
function stopcont_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stopcont (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object deletion, before destroying properties.
function stopcont_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to stopcont (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


