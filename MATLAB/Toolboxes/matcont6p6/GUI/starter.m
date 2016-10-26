function varargout = starter(varargin)
% STARTER Application M-file for starter.fig
%    FIG = STARTER launch starter GUI.
%    STARTER('callback_name', ...) invoke the named callback.

global gds path_sys MC driver_window FigPos;

if nargin == 0  % LAUNCH GUI
    file = fullfile(path_sys,gds.system,gds.diagram,gds.curve.new);
    f = strcat(file,'.mat');
    if exist(f,'file'),load(f);end    
    if ~exist('cds','var')    
        if isempty(gds.options), gds.options=contset;
            gds.options.PRC = 0;
            gds.options.dPRC = 0;
            gds.open.PRC = 0;
            gds.open.dPRC = 0;   
            MC.PRC = [];
            MC.dPRC = [];     
            gds.options.Input = 0;
        end
        gds.options=contset(gds.options,'Singularities',1);
        gds.options.IgnoreSingularity=[];
    end
    if ~isfield(gds.options,'PRC')
        gds.options.PRC = 0;
        gds.options.Input = 0;
    end
    if ~isfield(gds.options,'dPRC')
        gds.options.dPRC = 0;
        gds.options.Input = 0;
    end
    if ~gds.options.PRC
        MC.PRC = [];
        gds.open.PRC = 0;
    end
    if ~gds.options.dPRC
        MC.dPRC = [];
        gds.open.dPRC = 0;
    end
    gds.options.Multipliers = 1; 
    gds.options.Eigenvalues = 1;
    % XXXX
    gds.options.Locators = [];
    file= fullfile(path_sys,gds.system);
    save(file,'gds')
    file = strcat(file,'.mat');
    [ndim,r] = size(gds.parameters);
    load(file); 
    close(MC.starter);
    set(MC.mainwindow.compute,'Enable','on');
    set(MC.mainwindow.window,'Enable','on');
	fig = openfig(mfilename,'reuse');
    user.pos = [90.0 3 50 18]; user.num = 0;
    set(fig,'UserData',user);
	% Generate a structure of handles to pass to callbacks, and store it. 
   	handles  = guihandles(fig);
    user.pos = get(handles.slider1,'Position');
    set(handles.slider1,'Userdata',user);
    feval(gds.gui_starter,handles);
    guidata(fig, handles);
    MC.starter = fig;
    MC.starter_handles = guihandles(fig);
    set(fig,'Position',FigPos.starter);
    movegui(fig,'onscreen');
 	if nargout > 0
		varargout{1} = fig;
    end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
	try 
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard 
	catch
		errordlg(lasterr,mfilename);
        if ishandle(driver_window),delete(driver_window);end
	end

end

%--------------------------------------------------------------------------
function j = start_time(handles,j)
%enters the names of the coordinates in the window
%it the value is known, it is entered
global gds;
color = [1 1 1];
if strcmp(gds.time{1,1},'')
    string = '';
else
    string = num2str(gds.time{1,2});
end
j = j-2;
user.num = 0;
pos = [2 j 18 1.8];
user.pos = pos;
stat = uicontrol(handles.figuur,'Style','text','HorizontalAlignment','left','String',gds.time{1,1},'BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12,'UserData',user);
set(stat,'Position',pos,'UserData',user);
pos = [20 j 25 1.8];
user.pos = pos;
edit = uicontrol(handles.figuur,'Style','edit','HorizontalAlignment','left','String',string,'Tag','edit0','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(edit,'Callback','Callback');
set(edit,'Position',pos,'UserData',user);
guidata(handles.figuur,handles);


%--------------------------------------------------------------------------
function j = start_coordinates(handles,j)
%enters the names of the coordinates in the window
%it the value is known, it is entered
global gds;
color = [1 1 1];
user.num = 0;
for i = 1:(gds.dim)
    if strcmp(gds.coordinates{i,1},'')
        string = '';
    else
        string = num2str(gds.coordinates{i,2},'%0.8g');
    end
    tag1 = strcat('text',num2str(i));
    tag2 = strcat('edit',num2str(i));
    edit = uicontrol(handles.figuur,'Style','edit','HorizontalAlignment','left','String',string,'Tag',tag2,'BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
    set(edit,'Callback','Callback');
    j = j-2;
    pos = [2 j 18 1.8];
    user.pos = pos;
    stat = uicontrol(handles.figuur,'Style','text','HorizontalAlignment','left','String',gds.coordinates{i,1},'Tag',tag1,'BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
    set(stat,'Position',pos,'UserData',user);
    pos = [20 j 25 1.8];user.pos=pos;
    set(edit,'Position',pos,'UserData',user);
end
guidata(handles.figuur,handles);

%---------------------------------------------------------------------------
function j = start_parameters(handles,j)
%enters the names of the parameters in the window
%it the value is known, it is entered
global gds;
color = [1 1 1];
ndim = size(gds.parameters,1);
for i = 1:ndim   
    if strcmp(gds.parameters{i,1},'')
        string = '';
    else
        string = num2str(gds.parameters{i,2},'%0.8g');
    end
    value = 0;
    if ~isempty(gds.options.ActiveParams)
         str = gds.options.ActiveParams;
         x = str(str==i);
         if ~isempty(x)
             value = 1;
        end
    end
    tag2 = strcat('edit',num2str(i+gds.dim)); 
    user.num = num2str(i);
    edit = uicontrol(handles.figuur,'Style','edit','HorizontalAlignment','left','String',string,'Tag',tag2,'BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
    rad  = uicontrol(handles.figuur,'Style','radiobutton','String',gds.parameters{i,1},'Tag','ActiveParams','Callback','set_option','Max',1,'Value',value,'BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
    set(edit,'Callback','Callback');
    j = j-2;
    pos = [20 j 25 1.8];user.pos=pos;
    set(edit,'Position',pos,'UserData',user);
    pos = [2 j 18 1.8]; user.pos=pos;
    set(rad,'Position',pos,'UserData',user);
end
guidata(handles.figuur,handles);

%---------------------------------------------------------------------------
function j = start_parameters_orbit(handles,j)
%enters the names of the parameters in the window
%it the value is known, it is entered
global gds;
color = [1 1 1];
ndim = size(gds.parameters,1);
for i = 1:ndim   
    if strcmp(gds.parameters{i,1},'')
        string = '';
    else
        string = num2str(gds.parameters{i,2},'%0.8g');
    end
    tag2 = strcat('edit',num2str(i+gds.dim)); 
    user.num = num2str(i);
    edit = uicontrol(handles.figuur,'Style','edit','HorizontalAlignment','left','String',string,'Tag',tag2,'BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
    rad  = uicontrol(handles.figuur,'Style','text','HorizontalAlignment','left','String',gds.parameters{i,1},'Tag','ActiveParams','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
    set(edit,'Callback','Callback');
    j = j-2;
    pos = [20 j 25 1.8];user.pos=pos;
    set(edit,'Position',pos,'UserData',user);
    pos = [2 j 18 1.8]; user.pos=pos;
    set(rad,'Position',pos,'UserData',user);
end
guidata(handles.figuur,handles);
%--------------------------------------------------------------------------
function j = start_period(handles,j)
global gds;
dim = size(gds.parameters,1);
color = [1 1 1];
j = j-2;
post = [2 j 18 1.8];
user.num = 0; user.pos = post;
pose = [20 j 25 1.8];
string = num2str(gds.period,'%0.8g');
tag  = strcat('edit',num2str(gds.dim+dim+1));
stat = uicontrol(handles.figuur,'Style','text','HorizontalAlignment','left','String','Period','Tag','text_Period','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',post,'UserData',user);
edit = uicontrol(handles.figuur,'Style','text','HorizontalAlignment','left','String',gds.period,'Tag',tag,'Backgroundcolor',color,'units','characters','fontname','FixedWidth','fontsize',12);
user.pos = pose; set(edit,'Position',pose,'UserData',user);
guidata(handles.figuur,handles);

%-----------------------------------------------------------
function j = start_jacobian(handles,j)
global gds;
color = [1 1 1];
j = j-2;
user.num = 0;
pos  = [1 j 45 1.8]; user.pos=pos;
stat = uicontrol(handles.figuur,'Style','text','String','Jacobian Data','Tag','jacobian','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',pos,'UserData',user);
j = j-2;
post = [2 j 18 1.8];
pose = [20 j 25 1.8];
if isempty(gds.options.Increment)
    gds.options.Increment = 1e-5;
end
user.pos = post;
stat = uicontrol(handles.figuur,'Style','text','HorizontalAlignment','left','String','increment','Tag','text_increment','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',post,'UserData',user);
edit = uicontrol(handles.figuur,'Style','edit','HorizontalAlignment','left','String',num2str(gds.options.Increment),'Tag','Increment','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
user.pos = pose;set(edit,'Position',pose,'Callback','set_option','UserData',user);
guidata(handles.figuur,handles);

%--------------------------------------------------------------------------
function j = start_discret(handles,j)
%enter fields for 'ntst(number of testintervals' and 'ncol(number of collocationpoints'
global gds;
color = [1 1 1];
j = j-2;
pos = [1 j 45 1.8];
user.num = 0; user.pos=pos;
stat = uicontrol(handles.figuur,'Style','text','String','Discretization Data','Tag','discretization','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',pos,'UserData',user);
j = j-2;
post = [2 j 18 1.8];user.pos=post;
pose = [20 j 25 1.8];
stat = uicontrol(handles.figuur,'Style','text','HorizontalAlignment','left','String','ntst','Tag','text_ntst','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',post,'UserData',user);
edit = uicontrol(handles.figuur,'Style','edit','HorizontalAlignment','left','String',gds.discretization.ntst,'Tag','ntst','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
user.pos = pose;set(edit,'Position',pose,'Callback','set_option','UserData',user);
j = j-2;
post = [2 j 18 1.8];user.pos=post;
pose = [20 j 25 1.8];
stat = uicontrol(handles.figuur,'Style','text','HorizontalAlignment','left','String','ncol','Tag','text_ncol','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',post,'UserData',user);
edit = uicontrol(handles.figuur,'Style','edit','HorizontalAlignment','left','String',gds.discretization.ncol,'Tag','ncol','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
user.pos=pose;set(edit,'Position',pose,'Callback','set_option','UserData',user);
guidata(handles.figuur,handles);

%---------------------------------------------------------------
function j = start_multipliers(handles,j)
%enter the possibillity to compute the multipliers
global gds
color = [1 1 1];
j = j-2;
pos = [1 j 45 1.8];
user.num = 0; user.pos=pos;
stat = uicontrol(handles.figuur,'Style','text','String','Calculate multipliers','Tag','singularities','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',pos,'UserData',user);
j = j-2;
post = [2 j 22 1.8];user.pos=post;
pose = [24 j 21 1.8];        
stat = uicontrol(handles.figuur,'Style','text','HorizontalAlignment','left','String','multipliers','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
edit = uicontrol(handles.figuur,'Style','popupmenu','String',{'yes','no'},'Tag','Multipliers','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',post,'UserData',user);user.pos=pose;
set(edit,'Position',pose,'Callback','set_option','UserData',user);
guidata(handles.figuur,handles);


%---------------------------------------------------------------
function j = start_eigenvalues(handles,j)
%enter the possibillity to compute the eigenvalues
global gds
color = [1 1 1];
j = j-2;
pos = [1 j 45 1.8];user.num=0;user.pos=pos;
stat = uicontrol(handles.figuur,'Style','text','String','Calculate eigenvalues','Tag','eigenvalue','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',pos,'UserData',user);
j = j-2;
post = [2 j 18 1.8];user.pos=post;
pose = [22 j 23 1.8];        
stat = uicontrol(handles.figuur,'Style','text','HorizontalAlignment','left','String','eigenvalues','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
edit = uicontrol(handles.figuur,'Style','popupmenu','String',{'yes','no'},'Tag','Eigenvalues','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',post,'UserData',user);user.pos=pose;
set(edit,'Position',pose,'Callback','set_option','UserData',user);
guidata(handles.figuur,handles);

%---------------------------------------------------------------
function j = start_userfunctions(handles,j)
%enter the possibillity to compute the eigenvalues
global gds
if isempty(char(gds.userfunction)),   return;   end
dimu = size(gds.userfunction,2);
color = [1 1 1];
j = j-2;
pos = [0 j 40 1.8]; user.num = 0; user.pos = pos;
stat = uicontrol(handles.figuur,'Style','text','String','Monitor Userfunctions','Tag','userfunction','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',pos,'UserData',user);
for k= 1:dimu
    if ~isempty(gds.userfunction{k})
        tag1 = strcat('user',num2str(k));           
        string = gds.options.UserfunctionsInfo(k).name; 
        val = gds.options.UserfunctionsInfo(k).state;
        j=j-2; user.num = num2str(k);
        post = [2 j 20 1.8];user.pos = post;
        pose = [22 j 16 1.8];
        stat = uicontrol(handles.figuur,'Style','text','HorizontalAlignment','left','String',string,'Tag',tag1,'BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
        edit = uicontrol(handles.figuur,'Style','popupmenu','String',{'yes','no'},'Value',2-val,'Tag','UserfunctionsInfo','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
        set(stat,'Position',post,'UserData',user); user.pos = pose;
        set(edit,'Position',pose,'Callback','set_option','UserData',user);
    end
end
guidata(handles.figuur,handles);

% --------------------------------------------------------------------
function varargout = slider1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.slider1.
list = get(handles.figuur,'children');
step = get(h,'Value');
max  = get(h,'Max');
min = get(h,'min');
if (step>(max)|step==min)
    return
end
for i=1:length(list)
    if (strcmp(get(list(i),'Tag'),'slider1')~=1)
        user=get(list(i),'UserData');
        user.pos(2) = user.pos(2)+(max-step);
        set(list(i),'Position',user.pos);
    end
end   


%--------------------------------------------------------------------
function varargout = in_start(handles)
% Stub for Callback of the uicontrol handles.slider1.
list = get(handles.figuur,'children');
pos = get(list(1),'Position');
siz = pos(4);
max  = get(handles.slider1,'Max');
for i = 1:length(list)
    if (strcmp(get(list(i),'Tag'),'slider1')~=1)
        user = get(list(i),'UserData');
        user.pos(2) = user.pos(2)-(max-18);
        set(list(i),'Position',user.pos,'UserData',user);
    end
end   

%-------------------------------------------------------------------
function varargout = figuur_ResizeFcn(h, eventdata, handles, varargin)
% Get the figure size and position
global MC FigPos
Userd = get(h,'UserData');
if isempty(Userd),return;end
Orig_Size   = Userd.pos;%get the position of the figure before resize
Figure_Size = get(h,'Position');% get the current position
max = get(MC.starter_handles.slider1,'max');
min = get(MC.starter_handles.slider1,'min');
% If the resized figure is smaller than the 
% original figure size then compensate
if (Figure_Size(3) < 50)
    % Do not allow the width to become smaller
    % If the width is too small then reset to origianl width
    set(h,'Position',[Figure_Size(1) Figure_Size(2) Orig_Size(3) Figure_Size(4)]);
    Figure_Size = get(h,'Position');
    Userd.pos = Figure_Size;
    set(h,'UserData',Userd);
end 
if (Figure_Size(4) == Orig_Size(4))& (Figure_Size(4)<max)
    return;
end

if (Figure_Size(4) > max)
    % slider not visible
    set(handles.slider1,'Visible','off');
    set(h,'Position',[Figure_Size(1) Figure_Size(2) Figure_Size(3) max]);  
    siz = max-Orig_Size(4);
    list = get(h,'children');
    for i = 1:length(list)
        if (strcmp(get(list(i),'Tag'),'slider1')~=1)
            user = get(list(i),'UserData');
            user.pos(2) = user.pos(2)+siz;
            set(list(i),'Position',user.pos,'UserData',user);
        end
    end   
    Userd.pos=[Figure_Size(1) Figure_Size(2) Figure_Size(3) max];
    Figure_Size = Userd.pos;
    set(h,'UserData',Userd); 
elseif (Figure_Size(4) == max) 
    if Figure_Size(4)~=Orig_Size(4)
        set(h,'Position',[Figure_Size(1) Figure_Size(2) Figure_Size(3) max]);  
        siz = max-Orig_Size(4);
        list = get(h,'children');
        for i = 1:length(list)
            if (strcmp(get(list(i),'Tag'),'slider1')~=1)
                user = get(list(i),'UserData');
                user.pos(2) = user.pos(2)+siz;
                set(list(i),'Position',user.pos,'UserData',user);
            end
        end   
    end
    set(handles.slider1,'Visible','off');
    Userd.pos=[Figure_Size(1) Figure_Size(2) Figure_Size(3) max];
    Figure_Size = Userd.pos;
    set(h,'UserData',Userd); 
end

if (Figure_Size(4) > Orig_Size(4) && Figure_Size(4)< max)||(Figure_Size(4) < Orig_Size(4) && Figure_Size(4)>= min)
    % Does allow the height to change
    siz = Figure_Size(4)-Orig_Size(4);
    list = get(h,'children');
    for i = 1:length(list)
        if (strcmp(get(list(i),'Tag'),'slider1')~=1)
            user = get(list(i),'UserData');
            user.pos(2) = user.pos(2)+siz;
            set(list(i),'Position',user.pos,'UserData',user);
        end
    end   
    Userd.pos = Figure_Size;
    set(h,'UserData',Userd); 
    set(handles.slider1,'Visible','on');
end


%-------------------------------------------------------------------
function varargout = CloseStarter(varargin)
global MC FigPos;
if ishandle(MC.starter)
    FigPos.starter = get(MC.starter,'Position');
    delete(MC.starter);
end
MC.starter=[];


% XXXX
%---------------------------------------------------------------
function j = start_prc(handles,j)
%enter the possibillity to compute the phase response curve (for limit
%cycles)
global gds MC
color = [1 1 1];
j = j-2;
pos = [1 j 45 1.8];
user.num = 0; user.pos=pos;
stat = uicontrol(handles.figuur,'Style','text','String','Calculate Phase Response Curve','Tag','phase response curve','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',pos,'UserData',user);

j = j-2;
post = [2 j 20 1.8];user.pos=post;
pose = [20 j 25 1.6];  
if isempty(gds.options.Input)
    gds.options.Input = 0;
end
val = gds.options.Input;
stat = uicontrol(handles.figuur,'Style','text','HorizontalAlignment','left','String','Input','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
edit = uicontrol(handles.figuur,'Style','edit','HorizontalAlignment','left','String',num2str(val),'Tag','Input','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',post,'UserData',user);user.pos=pose;
set(edit,'Position',pose,'Callback','set_option','UserData',user);
guidata(handles.figuur,handles);

j = j-2;
post = [2 j 22 1.8];user.pos=post;
pose = [24 j 21 1.8]; 
if isempty(gds.options.PRC)
    gds.options.PRC = 0;
    MC.PRC = [];
    gds.open.PRC = 0;
end
val = gds.options.PRC;
if val > 0 && isempty(MC.PRC)
    MC.PRC = 1;
end
stat = uicontrol(handles.figuur,'Style','text','HorizontalAlignment','left','String','PRC','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
edit = uicontrol(handles.figuur,'Style','popupmenu','String',{'no','yes'},'Value',1+val,'Tag','PRC','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',post,'UserData',user);user.pos=pose;
set(edit,'Position',pose,'Callback','set_option','UserData',user);

j = j-2;
post = [2 j 22 1.8];user.pos=post;
pose = [24 j 21 1.8];   
if isempty(gds.options.dPRC)
    gds.options.dPRC = 0;
    MC.dPRC = [];
    gds.open.dPRC = 0;
end
val = gds.options.dPRC;
if val > 0 && isempty(MC.dPRC)
    MC.dPRC = 2;
end
stat = uicontrol(handles.figuur,'Style','text','HorizontalAlignment','left','String','dPRC','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
edit = uicontrol(handles.figuur,'Style','popupmenu','String',{'no','yes'},'Value',1+val,'Tag','dPRC','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',post,'UserData',user);user.pos=pose;
set(edit,'Position',pose,'Callback','set_option','UserData',user);

