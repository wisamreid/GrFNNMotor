function varargout = numeric(varargin)
% NUMERIC Application M-file for numeric.fig
%    FIG = NUMERIC launch numeric GUI.
%    NUMERIC('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 14-Nov-2001 12:21:37
global gds path_sys MC FigPos cds
file = fullfile(path_sys,gds.system);save(file,'gds');
if nargin ==0  % LAUNCH GUI
    file = strcat(file,'.mat');
    load(file);
    close(MC.numeric_fig);
    fig = openfig(mfilename,'reuse');
    user=[];
    user.num=0;
    user.pos=[120 3. 50.0000  15];
	% Generate a structure of handles to pass to callbacks, and store it. 
	handles  = guihandles(fig);
    set(fig,'UserData',user);   
    
    user.pos = get(handles.slider1,'Position');
    set(handles.slider1,'Userdata',user);
    feval(gds.gui_numeric,handles); 
    guidata(fig,handles);
    MC.numeric_fig = fig;
    MC.numeric_handles = guihandles(fig);
    set(fig,'Position',user.pos);
    set(fig,'Position',FigPos.numeric_fig);
    
    if verLessThan('matlab' , '8.4.0')
        %do nothing
    else
        set(fig, 'ResizeFcn', @(h,e) 0);
        set(fig, 'SizeChangedFcn', @(h,e)  resizedelay( @() numeric('figuur_ResizeFcn',h,e,guidata(h))));
    end
    
    movegui(fig,'onscreen');
  	if nargout > 0
		varargout{1} = fig;
	end
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
        guidata(handles.numeric_fig,handles);
    catch
	end

end

%--------------------------------------------------------------------------
function j = start_time(handles,j)
%enters the names of the coordinates in the window
%it the value is known, it is entered
global gds 
color = [1 1 1];
stat1 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String',gds.time{1,1},'BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'units','characters');
stat2 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','','Tag','time','BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'units','characters');
j = j-2;
pos = [2 j 18 1.8];user.num=0;user.pos=pos;
set(stat1,'Position',pos,'UserData',user);
pos = [20 j 27 1.8];user.pos=pos;
set(stat2,'Position',pos,'UserData',user);
% guidata(handles.numeric_fig,handles);
% MC.numeric_handles=guihandles(fig)
%--------------------------------------------------------------------------
function j=start_coordinates(handles,j)
global gds 
color = [1 1 1];
j = j-2;user.num=0;
stat = uicontrol(handles.numeric_fig,'Style','text','String','Coordinates','BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
pos = [12 j 18 1.8];user.pos=pos;
set(stat,'Position',pos,'UserData',user);
for i = 1:(gds.dim)
    string = '';
    tag1 = strcat('text',num2str(i));
    tag2 = strcat('coord',num2str(i));    
    stat1 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String',gds.coordinates{i,1},'Tag',tag1,'BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
    stat2 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String',string,'Tag',tag2,'BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
    j = j-2;
    pos = [2 j 18 1.8];user.pos=pos;
    set(stat1,'Position',pos,'UserData',user);    
    pos = [20 j 27 1.8];user.pos=pos;
    set(stat2,'Position',pos,'UserData',user);
end
% guidata(handles.numeric_fig,handles);

%---------------------------------------------------------------------------
function j = start_parameters(handles,j)
global gds
color = [1 1 1];
j = j-2;
stat = uicontrol(handles.numeric_fig,'Style','text','String','Parameters','BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
pos = [12 j 18 1.8];user.num=0;user.pos=pos;
set(stat,'Position',pos,'UserData',user);
ndim = size(gds.parameters,1);
for i=1:ndim
    string='';
    tag1 = strcat('text2',num2str(i));
    tag2 = strcat('param',num2str(i));
    stat1 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String',gds.parameters{i,1},'Tag',tag1,'BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
    stat2 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String',string,'Tag',tag2,'BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
    j = j-2;
    poss = [2 j 18 1.8];user.pos=poss;
    set(stat1,'Position',poss,'Userdata',user);
    pos = [20 j 27 1.8];user.pos=pos;
    set(stat2,'Position',pos,'UserData',user);
end
% guidata(handles.numeric_fig,handles);

%--------------------------------------------------------------------------
function j = start_period(handles,j)
global gds MC ;
dim = size(gds.parameters,1);
color = [1 1 1];
j = j-2;
stat = uicontrol(handles.numeric_fig,'Style','text','String','Period','BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
pos = [12 j 18 1.8];user.num=0;user.pos=pos;
set(stat,'Position',pos,'UserData',user);
j = j-2;
post = [2 j 18 1.8];user.pos=post;
pose = [20 j 27 1.8];
stat1 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','Period','Tag','text_Period','BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
set(stat1,'Position',post,'UserData',user);user.pos=pose;
stat2 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','','Tag','period','Backgroundcolor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
set(stat2,'Position',pose,'UserData',user);
% guidata(handles.numeric_fig,handles);

%-----------------------------------------------------------
function j = start_multipliers(handles,j)
global gds 
color = [1 1 1];
j = j-2;
pos = [12 j 18 1.8];user.num=0;user.pos=pos;
stat = uicontrol(handles.numeric_fig,'Style','text','String','Multipliers','Tag','multipliers','BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
set(stat,'Position',pos,'UserData',user);
jarg = j -2*gds.dim;
for k = 1:gds.dim
    j = j-2;  
    tag1 = sprintf('Mod[%d]',k);
    tag11=strcat('Mod_',num2str(k));
    tag2 = sprintf('Arg[%d]',k);
    tag22 = strcat('Arg_',num2str(k));
    post = [1 j 18 1.8];user.pos=post;
    pose = [19 j 27 1.8];             
    stat1 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String',tag1,'BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
    stat2 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','','Tag',tag11,'BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
    set(stat1,'Position',post,'UserData',user);user.pos=pose;
    set(stat2,'Position',pose,'UserData',user);
    jarg = jarg-2;
    post = [1 jarg 18 1.8];user.pos=post;
    pose = [19 jarg 27 1.8];             
    stat1 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String',tag2,'BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
    stat2 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','','Tag',tag22,'BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
    set(stat1,'Position',post,'UserData',user);user.pos=pose;
    set(stat2,'Position',pose,'UserData',user);
end
j=jarg;
% guidata(handles.numeric_fig,handles);

%---------------------------------------------------------------
function j = start_userfunctions(handles,j)
%enter the possibillity to compute the eigenvalues
global gds 
color = [1 1 1];
j = j-2;
pos = [10 j 28 1.8]; user.num = 0; user.pos = pos;
stat = uicontrol(handles.numeric_fig,'Style','text','String','User functions','Tag','userfunction','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',pos,'UserData',user);
if isfield(gds,'userfunction')
    dimu = size(gds.userfunction,2);
else return
end
for k= 1:dimu
    if (gds.options.UserfunctionsInfo(k).state==0)
        continue;
    end
    tag1 = strcat('user_',num2str(k));           
    string = gds.options.UserfunctionsInfo(k).name; 
    j = j-2; user.num = num2str(k);
    post = [2 j 20 1.8];user.pos = post;
    pose = [22 j 17 1.8];
    stat1 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String',string,'Tag','user','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
    stat2 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','','Tag',tag1,'BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
    set(stat1,'Position',post,'UserData',user); user.pos = pose;
    set(stat2,'Position',pose,'UserData',user);
end
% guidata(handles.numeric_fig,handles);

%---------------------------------------------------------------
function j = start_stepsize(handles,j)
global gds MC;
color = [1 1 1];
j = j-2;
post = [2 j 18 1.8];user.num=0;user.pos=post;
pose = [20 j 27 1.8];             
stat1 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','Stepsize','BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
stat2 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','','Tag','stepsize','BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
set(stat1,'Position',post,'UserData',user);user.pos=pose;
set(stat2,'Position',pose,'UserData',user);
% guidata(handles.numeric_fig,handles);

%---------------------------------------------------------------
function j = start_npoints(handles,j)
global gds 
color = [1 1 1];
j = j-2;
post = [2 j 18 1.8];user.num=0;user.pos=post;
pose = [20 j 22 1.8];             
stat1 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','Npoints','BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
stat2 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','','Tag','npoints','BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
set(stat1,'Position',post,'UserData',user);user.pos=pose;
set(stat2,'Position',pose,'UserData',user);
% guidata(handles.numeric_fig,handles);
%-----------------------------------------------------------
function j = start_eigenvalues(handles,j)
global gds MC 
color = [1 1 1];
j = j-2;
pos = [12 j 18 1.8];user.num=0;user.pos=pos;
stat = uicontrol(handles.numeric_fig,'Style','text','String','Eigenvalues','Tag','eigenvalues','BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
set(stat,'Position',pos,'UserData',user);
jim = j -2*gds.dim;
for k = 1:gds.dim
    j = j-2;    
    tag1 = sprintf('Re[%d]',k);
    tag11 = sprintf('Re_%d',k);
    tag2 = sprintf('Im[%d]',k);
    tag22 = sprintf('Im_%d',k);
    post = [1 j 18 1.8];user.pos=post;
    pose = [19 j 27 1.8];             
    stat1 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String',tag1,'BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
    stat2 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','','Tag',tag11,'BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
    set(stat1,'Position',post,'UserData',user);user.pos=pose;
    set(stat2,'Position',pose,'UserData',user);
    jim=jim-2;
    post = [1 jim 18 1.8];user.pos=post;
    pose = [19 jim 27 1.8];             
    stat1 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String',tag2,'BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
    stat2 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','','Tag',tag22,'BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
    set(stat1,'Position',post,'UserData',user);user.pos=pose;
    set(stat2,'Position',pose,'UserData',user);

end
j=jim;
% guidata(handles.numeric_fig,handles);
% MC.numeric_handles = handles;

% --------------------------------------------------------------------
function slider1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.slider1.
list = get(handles.numeric_fig,'children');
step = get(h,'Value');
max  = get(h,'Max');min=get(h,'min');
if (step > max||step==min)
    return
end
for i = 1:length(list)
    if (strcmp(get(list(i),'Tag'),'slider1')~=1)&&(strcmp(get(list(i),'Tag'),'window')~=1)
        user = get(list(i),'UserData');
        user.pos(2) = user.pos(2)+(max-step);
        set(list(i),'Position',user.pos);
    end
end   

% --------------------------------------------------------------------
function varargout = in_numeric(handles)
% Stub for Callback of the uicontrol handles.slider1.
list = get(handles.numeric_fig,'children');
max  = get(handles.slider1,'Max');
for i=1:length(list)
    if (strcmp(get(list(i),'Tag'),'slider1')~=1)&&(strcmp(get(list(i),'Tag'),'window')~=1)
        user = get(list(i),'UserData');
        user.pos(2) = user.pos(2)-(max-15);
        set(list(i),'Position',user.pos,'UserData',user);
    end
end   
%-------------------------------------------------------------------
function figuur_ResizeFcn(h, eventdata, handles, varargin)
% Get the figure size and position
Userd = get(h,'UserData');
if isempty(Userd),return;end
Orig_Size   = Userd.pos;%get the position of the figure before resize
Figure_Size = get(h,'Position');% get the current position
max = get(handles.slider1,'max');
min = get(handles.slider1,'min');
% If the resized figure is smaller than the 
% original figure size then compensate

if verLessThan('matlab' , '8.4.0') %NN fix 2014b
    if (Figure_Size(3) ~= Orig_Size(3))
        % Do not allow the width to change
        set(h,'Position',[Figure_Size(1) Figure_Size(2) Orig_Size(3) Figure_Size(4)]);
        Figure_Size(3) = Orig_Size(3);
        Userd.pos=Figure_Size;
        set(h,'UserData',Userd);
    end
end

if (Figure_Size(4) == Orig_Size(4))&& (Figure_Size(4)<max)
    return;
end

if (Figure_Size(4)== max)
    if Figure_Size(4)~=Orig_Size(4)
        siz = max-Orig_Size(4); 
        list = get(h,'children');
        for i = 1:length(list)
            if (strcmp(get(list(i),'Tag'),'slider1')==1)
                user = get(list(i),'UserData');
                user.pos(4) = user.pos(4)+siz;
                set(list(i),'Position',user.pos,'UserData',user);
            elseif (strcmp(get(list(i),'Tag'),'window')~=1)
                user = get(list(i),'UserData');
                user.pos(2) = user.pos(2)+siz;
                set(list(i),'Position',user.pos,'UserData',user);
            end
        end       
    end
    %Do not allow the height to change
    set(handles.slider1,'Visible','off');
    Userd.pos=[Figure_Size(1) Figure_Size(2) Figure_Size(3) max];
    Figure_Size = Userd.pos;
    set(h,'UserData',Userd); 
    set(h,'Position',Figure_Size);  
    return
elseif (Figure_Size(4) >max)
    % Do not allow the height to change
    set(handles.slider1,'Visible','off');
    Userd.pos=[Figure_Size(1) Figure_Size(2) Figure_Size(3) max];
    Figure_Size = Userd.pos;
    set(h,'UserData',Userd); 
    set(h,'Position',Figure_Size);  
    siz = max-Orig_Size(4);
    list = get(h,'children');
    for i = 1:length(list)
         if (strcmp(get(list(i),'Tag'),'slider1')==1)
             user = get(list(i),'UserData');
             user.pos(4) = user.pos(4)+siz;
             set(list(i),'Position',user.pos,'UserData',user);
        elseif (strcmp(get(list(i),'Tag'),'window')~=1)
            user = get(list(i),'UserData');
            user.pos(2) = user.pos(2)+siz;
            set(list(i),'Position',user.pos,'UserData',user);
        end
    end     
    return
end
if (Figure_Size(4) > Orig_Size(4) && Figure_Size(4)< max)||(Figure_Size(4) < Orig_Size(4) && Figure_Size(4)>= min)
    % Does allow the height to change
    siz = Figure_Size(4)-Orig_Size(4);
    list = get(h,'children');
    for i = 1:length(list)
        if (strcmp(get(list(i),'Tag'),'slider1')==1)
            user = get(list(i),'UserData');
            user.pos(4) = user.pos(4)+siz;
            set(list(i),'Position',user.pos,'UserData',user);
        elseif (strcmp(get(list(i),'Tag'),'window')~=1)
            user = get(list(i),'UserData');
            user.pos(2) = user.pos(2)+siz;
            set(list(i),'Position',user.pos,'UserData',user);
        end
    end   
    Userd.pos = Figure_Size;
    set(h,'UserData',Userd); 
    set(handles.slider1,'Visible','on');
end

if ((Figure_Size(4)<min) && verLessThan('matlab' , '8.4.0')) %fix 2014b (NN)
    % Do not allow the height to change
    set(h,'Position',Orig_Size); 
end


%----------------------------------------------------------------------
function CloseNumeric(varargin)
global MC calculation_progress FigPos;%lds;
if calculation_progress~=0 
   return;
end;
if ishandle(MC.numeric_fig)
    FigPos.numeric_fig = get(MC.numeric_fig,'Position');
    delete(MC.numeric_fig);
end
MC.numeric_fig = [];
MC.numeric_handles=[];
