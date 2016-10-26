function varargout = continuer(varargin)
% CONTINUER Application M-file for continuer.fig
%    FIG = CONTINUER launch continuer GUI.
%    CONTINUER('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 04-Mar-2002 14:52:18
global gds MC driver_window FigPos;
if nargin == 0  % LAUNCH GUI
    close(MC.continuer);
    close(MC.integrator);
	fig = openfig(mfilename,'reuse');
    set(fig,'Position',FigPos.continuer);
    movegui(fig,'onscreen');
   	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);  
    MC.continuer = fig;
    if ~isempty(gds.options.MinStepsize)
        set(handles.MinStepsize,'String',gds.options.MinStepsize);
    end
    if ~isempty(gds.options.MaxStepsize)
        set(handles.MaxStepsize,'String',gds.options.MaxStepsize);
    end
    if ~isempty(gds.options.InitStepsize)
        set(handles.InitStepsize,'String',gds.options.InitStepsize);
    end
    if ~isempty(gds.options.FunTolerance)
        set(handles.FunTolerance,'String',gds.options.FunTolerance);
    end
    if ~isempty(gds.options.VarTolerance)
         set(handles.VarTolerance,'String',gds.options.VarTolerance);
    end
    if ~isempty(gds.options.TestTolerance)
        set(handles.TestTolerance,'String',gds.options.TestTolerance);
    end
    if ~isempty(gds.options.MaxNewtonIters)
        set(handles.MaxNewtonIters,'String',gds.options.MaxNewtonIters);
    end
    if ~isempty(gds.options.MaxCorrIters)
        set(handles.MaxCorrIters,'String',gds.options.MaxCorrIters);
    end
    if ~isempty(gds.options.MaxTestIters)
        set(handles.MaxTestIters,'String',gds.options.MaxTestIters);
    end
    if ~isempty(gds.options.MaxNumPoints)
        set(handles.MaxNumPoints,'String',gds.options.MaxNumPoints);
    end
    if ~isempty(gds.options.CheckClosed)
        set(handles.CheckClosed,'String',gds.options.CheckClosed);
    end    
    if ~isempty(gds.options.Adapt)
        set(handles.Adapt,'String',gds.options.Adapt);
    end
    guidata(fig, handles);

	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
		errordlg(lasterr,'error continuer');
        if ishandle(driver_window),delete(driver_window);end
	end

end

%-----------------------------------------------------------------------
function CloseContinuer(varargin)
global MC FigPos;
    if ishandle(MC.continuer)
        FigPos.continuer = get(MC.continuer,'Position');
        delete(MC.continuer);
    end
    MC.continuer=[];

%-------------------------------------------------------------------
function figuur_ResizeFcn(h, eventdata, handles, varargin)
% Get the figure size and position
global MC FigPos
Userd = get(h,'UserData');
if isempty(Userd)
    Orig_Size = get(h,'Position');
else
    Orig_Size   = Userd.pos;%get the position of the figure before resize
end
Figure_Size = get(h,'Position');% get the current position
% % If the resized figure is smaller than the 
% % original figure size then compensate
if (Figure_Size(3) < 50)
    % Do not allow the width to become smaller
    % If the width is too small then reset to origianl width
    set(h,'Position',[Figure_Size(1) Figure_Size(2) 50 Figure_Size(4)]);
    Figure_Size = get(h,'Position');
else
    set(h,'Position',[Figure_Size(1) Figure_Size(2) Figure_Size(3) Figure_Size(4)]);
    Figure_Size = get(h,'Position');
end
Userd.pos = Figure_Size;
set(h,'UserData',Userd);
