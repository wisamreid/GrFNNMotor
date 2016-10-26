function varargout = integrator(varargin)
% INTEGRATOR Application M-file for integrator.fig
%    FIG = INTEGRATOR launch integrator GUI.
%    INTEGRATOR('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 18-Feb-2002 16:46:12
global gds path_sys MC driver_window FigPos;
if nargin == 0  % LAUNCH GUI
    close(MC.continuer);
    close(MC.integrator);
    set(MC.mainwindow.compute,'Enable','on');
    fig = openfig(mfilename,'reuse');
    set(fig,'Position',FigPos.integrator);
    movegui(fig,'onscreen');
	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
    MC.integrator = fig;
    load_integrator(handles);

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
function method_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.method.
global gds path_sys;
index=get(handles.method,'Value');
switch index
case 1
    gds.integrator.method='ode45';
    set(handles.Refine,'String','4');
    if isfield(handles,'BDF') && ishandle(handles.BDF)
        delete(handles.BDF);
        delete(handles.MaxOrder);
        delete(handles.BDF1);
        delete(handles.MaxOrder1);
    end   
case 2
    gds.integrator.method='ode23';
    set(handles.Refine,'String','1');
    if isfield(handles,'BDF') && ishandle(handles.BDF)      
        delete(handles.BDF);
        delete(handles.MaxOrder);
        delete(handles.BDF1);
        delete(handles.MaxOrder1);
    end
case 3
    gds.integrator.method='ode113';  
    set(handles.Refine,'String','1');
    if isfield(handles,'BDF') && ishandle(handles.BDF)
        delete(handles.BDF);
        delete(handles.MaxOrder);
        delete(handles.BDF1);
        delete(handles.MaxOrder1);
    end  
case 4
    gds.integrator.method='ode15s';
    stat=uicontrol(handles.integrator,'Style','text','String','BDF','tag','BDF1','HorizontalAlignment','left','BackGroundColor',[1 1 1],'units','normalized','fontname','FixedWidth','fontsize',12);
    edit=uicontrol(handles.integrator,'Style','popupmenu','String',{'Off','On'},'Tag','BDF','BackGroundColor',[1 1 1],'units','normalized','fontname','FixedWidth','fontsize',12);
    if isfield(gds.integrator.options,'BDF')&&~isempty(gds.integrator.options.BDF)
        str=gds.integrator.options.BDF;
        switch str
        case 'On',set(edit,'Value',2);
        otherwise set(edit,'Value',1);            
        end
    end
    set(edit,'Callback','integrator(''BDF_callback'',guidata(gcbo))');
    pos=[0.023 0.135 0.55 0.075];
    set(stat,'Position',pos);
    pos=[0.55 0.135 0.40 0.075];
    set(edit,'Position',pos);
    stat=uicontrol(handles.integrator,'Style','text','String','MaxOrder','tag','MaxOrder1','HorizontalAlignment','left','BackGroundColor',[1 1 1],'units','normalized','fontname','FixedWidth','fontsize',12);
    edit=uicontrol(handles.integrator,'Style','popupmenu','String',{'5','4','3','2','1'} ,'tag','MaxOrder','BackGroundColor',[1 1 1],'units','normalized','fontname','FixedWidth','fontsize',12); 
    if isfield(gds.integrator.options,'MaxOrder')&&~isempty(gds.integrator.options.MaxOrder)
        str=gds.integrator.options.MaxOrder;
        switch str
        case '1',set(edit,'Value',5);
        case '4',set(edit,'Value',2);
        case '3',set(edit,'Value',3);
        case '2',set(edit,'Value',4);
        otherwise set(edit,'Value',1);            
        end
    end
    set(edit,'Callback','integrator(''MaxOrder_callback'',guidata(gcbo))');
    pos=[0.023 0.048 0.55 0.075];
    set(stat,'Position',pos);
    pos=[0.55 0.048 0.40 0.075];
    set(edit,'Position',pos);
    set(handles.Refine,'String','1');
    handles=guihandles(handles.integrator);
    guidata(handles.integrator,handles);
case 5
    gds.integrator.method='ode23s';
    set(handles.Refine,'String','1');
     if isfield(handles,'BDF') && ishandle(handles.BDF)
        delete(handles.BDF);
        delete(handles.MaxOrder);
        delete(handles.BDF1);
        delete(handles.MaxOrder1);
    end 
case 6
    gds.integrator.method='ode23t';
    set(handles.Refine,'String','1');
    if isfield(handles,'BDF') && ishandle(handles.BDF)
        delete(handles.BDF);
        delete(handles.MaxOrder);
        delete(handles.BDF1);
        delete(handles.MaxOrder1);
    end 
case 7
    gds.integrator.method='ode23tb';
    set(handles.Refine,'String','1');
    if isfield(handles,'BDF') && ishandle(handles.BDF)
        delete(handles.BDF);
        delete(handles.MaxOrder);
        delete(handles.BDF1);
        delete(handles.MaxOrder1);
    end 
case 8
    gds.integrator.method='ode78';
    set(handles.Refine,'String','1');
    if isfield(handles,'BDF') && ishandle(handles.BDF)
        delete(handles.BDF);
        delete(handles.MaxOrder);
        delete(handles.BDF1);
        delete(handles.MaxOrder1);
    end   
case 9
    gds.integrator.method='ode87';
    set(handles.Refine,'String','1');
    if isfield(handles,'BDF') && ishandle(handles.BDF)
        delete(handles.BDF);
        delete(handles.MaxOrder);
        delete(handles.BDF1);
        delete(handles.MaxOrder1);
    end   
end
file=fullfile(path_sys,gds.system);
save(file,'gds');

% --------------------------------------------------------------------
function interval_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.interval.
global gds path_sys;
val=str2double(get(handles.interval,'String'));
if val<0 then
    errordlg('this isn''t possible');
    set(handles.interval,'string','1');
else
    gds.integrator.tspan=[gds.time{1,2} (val+gds.time{1,2})];
end
file=fullfile(path_sys,gds.system);
save(file,'gds');


% --------------------------------------------------------------------
function MaxStepsize_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.MaxStepsize.
global gds path_sys;
val=str2double(get(handles.MaxStepsize,'String'));
gds.integrator.options=odeset(gds.integrator.options,'MaxStep',val);
file=fullfile(path_sys,gds.system);
save(file,'gds');


% --------------------------------------------------------------------
function InitStepsize_Callback(h, eventdata, handles, varargin)
global gds path_sys;
val=str2double(get(handles.InitStepsize,'String'));
gds.integrator.options=odeset(gds.integrator.options,'InitialStep',val);
file=fullfile(path_sys,gds.system);
save(file,'gds');


% --------------------------------------------------------------------
function RelTolerance_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.RelTolerance.
global gds path_sys;
val=str2double(get(handles.RelTolerance,'String'));
gds.integrator.options=odeset(gds.integrator.options,'RelTol',val);
file=fullfile(path_sys,gds.system);
save(file,'gds');


% --------------------------------------------------------------------
function edit5_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.edit5.
global gds path_sys;
val=str2double(get(handles.edit5,'String'));
gds.integrator.options=odeset(gds.integrator.options,'AbsTol',val);
file=fullfile(path_sys,gds.system);
save(file,'gds');


% --------------------------------------------------------------------
function Refine_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.Refine.
global gds path_sys;
val=str2double(get(handles.Refine,'String'));
gds.integrator.options=odeset(gds.integrator.options,'Refine',val);
file=fullfile(path_sys,gds.system);
save(file,'gds');


% --------------------------------------------------------------------
function Normcontrol_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.Normcontrol.
global gds path_sys;
gds.integrator.options=[];
val=get(handles.Normcontrol,'value');
switch val
case 1
    str='off';
case 2
    str='on';
end
gds.integrator.options=odeset(gds.integrator.options,'NormControl',str);
file=fullfile(path_sys,gds.system);
save(file,'gds');


%-------------------------------------------------------------------
function BDF_callback(handles)
global gds path_sys;
index=get(handles.BDF,'Value');
switch index
case 1
    gds.integrator.options=odeset(gds.integrator.options,'BDF','off');
case 2
    gds.integrator.options=odeset(gds.integrator.options,'BDF','on');
end
file=fullfile(path_sys,gds.system);
save(file,'gds');

%-----------------------------------------------------------------
function MaxOrder_callback(handles)
global gds path_sys;
index=get(handles.MaxOrder,'Value');
string=get(handles.MaxOrder,'String');
gds.integrator.options=odeset(gds.integrator.options,'MaxOrder',str2double(string{index}));
file=fullfile(path_sys,gds.system);
save(file,'gds');


%-----------------------------------------------------------------
function load_integrator(handles)
global gds;
%not DO possible for ode78,ode87 
if strcmp(gds.type,'DO ')
    str_method=get(handles.method,'String');
    str_method(8:9)=[];
    set(handles.method,'String',str_method);
end
switch gds.integrator.method
case 'ode45'
    val=1;
case 'ode23'
    val=2;
case 'ode113';
    val=3;
case 'ode15s'
    val=4;
    set(handles.method,'Value',val);
    integrator('method_Callback',handles.integrator,[],guidata(handles.integrator));    
case 'ode23s'
    val=5;
case 'ode23t'
    val=6;
case 'ode78'
    val=8;
case 'ode87'
    val=9;
otherwise
    val=7;
end
set(handles.method,'Value',val);
set(handles.interval,'string',num2str(abs(gds.integrator.tspan(2)-gds.integrator.tspan(1))));
if isfield(gds.integrator.options,'InitialStep')&&~isempty(gds.integrator.options.InitialStep)
    str=gds.integrator.options.InitialStep;
else
    str='<automatic>';
end
set(handles.InitStepsize,'string',str);
if isfield(gds.integrator.options,'MaxStep')&&~isempty(gds.integrator.options.MaxStep)
    str=gds.integrator.options.MaxStep;
else
    str='<automatic>';
end
set(handles.MaxStepsize,'string',str);
if isfield(gds.integrator.options,'RelTol')&&~isempty(gds.integrator.options.RelTol)
    str=gds.integrator.options.RelTol;
else
    str='1e-3';
end
set(handles.RelTolerance,'string',str);
if isfield(gds.integrator.options,'AbsTol')&&~isempty(gds.integrator.options.AbsTol)
    str=gds.integrator.options.AbsTol;
else
    str='1e-6';
end
set(handles.edit5,'string',str);

if isfield(gds.integrator.options,'Refine')&&~isempty(gds.integrator.options.Refine)
    str=gds.integrator.options.AbsTol;
else
    str='1';
end
set(handles.Refine,'string',str);
if isfield(gds.integrator.options,'Refine')&&~isempty(gds.integrator.options.Refine)
    str=gds.integrator.options.AbsTol;
else
    str='1';
end
set(handles.Refine,'string',str);
if isfield(gds.integrator.options,'NormControl')&&~isempty(gds.integrator.options.NormControl)
    str=gds.integrator.options.NormControl;
    switch str
        case 'on'
            set(handles.Normcontrol,'Value',2);
        otherwise 
            set(handles.Normcontrol,'Value',1);
    end
end
%------------------------------------------------------------------

function CloseIntegrator(varargin)
global MC FigPos;
    if ishandle(MC.integrator)
        FigPos.integrator = get(MC.integrator,'Position');
        delete(MC.integrator);
    end
    MC.integrator = [];
