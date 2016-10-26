function varargout = select_diagram(varargin)
% EDITSYSTEM Application M-file for editsystem.fig
%    FIG = EDITSYSTEM launch editsystem GUI.
%    EDITSYSTEM('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 20-Sep-2002 13:41:09
global gds path_sys driver_window
if nargin <= 0  % LAUNCH GUI
    if nargin==0
        initial_dir=fullfile(path_sys,gds.system);
    elseif nargin==1 && exist(varargin{1},'dir')
        initial_dir=varargin{1};
    else 
        errordlg('Input argument must be a valid directory',...
            'Input Argument Error!')
        return
    end
    
    fig = openfig(mfilename,'reuse');

	% Use system color scheme for figure:
	set(fig,'Color',get(0,'defaultUicontrolBackgroundColor'));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
    %Populate the listbox
    load_listbox(initial_dir,handles);
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
function select_diagram_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.files.
global gds path_sys MC FigPos;
initial_dir = fullfile(path_sys,gds.system);
index_selected = get(handles.listbox_diagram,'Value');
dir_list = get(handles.listbox_diagram,'String');
dia = dir_list{index_selected};
diagramname = fullfile(initial_dir,dia);
cd (diagramname);w = what;
gds.diagram = dia;
if size(w.mat,1)>1
    filename = fullfile(initial_dir,dia,w.mat{1});
    load(filename); index = 1;
    if exist('cds','var')
        cds.options.UserfunctionsInfo = gds.options.UserfunctionsInfo;
        cds.options.Userfunctions = gds.options.Userfunctions;
        gds.options = cds.options;
    end
    gds.type = ctype; gds.point = point;feval(strcat('gui_',deblank(ctype)));
    
            sn = w.mat{1};
            if size(sn,1) > 1
                sn = sn(1,:);
            end
    gds.curve.new = strrep(sn,'.mat',''); gds.curve.old = strrep(sn,'.mat','');
    file = fullfile(path_sys,gds.system); save(file,'gds');
    feval(gds.gui_load_point,index,x,'not_save',filename,1);
    list=get(MC.mainwindow.initial_point,'children');
    for i=1:length(list)
        tag=get(list(i),'Tag');[tag1,tag2]=strtok(tag,'_');
        if strcmp(deblank(point),tag1)
            break
        end
    end
    [point,str] = strtok(tag,'_');%point=EP;str=EP
    set(MC.mainwindow.window,'enable','on');
    list = get(MC.mainwindow.curvetype,'children');%find all the type of curves
    for i = 1:length(list)
        type = get(list(i),'tag');
        if isempty(findstr(str,type))
            set(list(i),'enable','off');%makes all other continuation of curves not possible
        else
            set(list(i),'enable','on');%EP-curve is possible
        end
    end  
    if exist('cds','var');continuer;else integrator;end
    starter;%load the appropriate continuer and starter window (see starter)
else
    gds.point = '';gds.type = '';gds.curve.new = ''; 
    set(MC.mainwindow.compute,'enable','off');
    set(MC.mainwindow.window,'enable','off');
    set(MC.mainwindow.Type,'enable','on');
    set(MC.mainwindow.select_userfunctions,'enable','on');
    %starter-window open?
    close(MC.starter); gds.open.figuur = 0;
    % continuer-window open?
    close(MC.continuer); gds.open.continuer = 0;
    % numeric window open?
    close(MC.numeric_fig);gds.open.numeric_fig = 0;
    %2D-plot open
    close(MC.D2); gds.open.D2 = 0;
    %3D-plot open
    close(MC.D3); gds.open.D3 = 0;
    %PRC-plot open
    close(MC.PRC); 
    MC.PRC = [];
    gds.open.PRC = 0;
    close(MC.integrator); gds.open.integrator = 0;
    %PRC-plot open
    close(MC.dPRC); 
    MC.dPRC = [];gds.open.dPRC = 0;
end
load_matcont;
close(handles.diagramfig);
cd(path_sys);cd ..;
    
%-----------------------------------------------------------------
function load_listbox(dir_path,handles)
global path_sys gds;
string='';val=1;
if ~(exist(dir_path,'dir')~=7)
    cd (dir_path);d = dir;
    ind = find(vertcat(d.isdir));
    if size(ind,1)>2; 
        string = cellstr(char(d(ind(3:end,1)).name));
    end
    val = strmatch(gds.diagram,string,'exact');
end
if isempty(val), val = 1;end
guidata(handles.diagramfig,handles);
set(handles.listbox_diagram,'String',string,'Value',val);
cd(path_sys);cd ..;

%-------------------------------------------------------------------
function actions_delete_Callback(h,eventdata,handles,varargin)
global gds path_sys MC;
index_selected = get(handles.listbox_diagram,'Value');
file_list = get(handles.listbox_diagram,'String');
dir1 = file_list{index_selected};
dir1 = fullfile(path_sys,gds.system,dir1);
cd(dir1);delete *.*;
cd(path_sys);cd ..;
dir=fullfile(path_sys,gds.system);
gds.diagram='diagram';
cd(dir);
str=sprintf('1)Select the directory named ''%s''to remove in the following ''delete box''\n2)Right-click and select Delete from the context menu.\n3)Close the ''delete box''',dir1);
button = questdlg(str,'Info deleting diagram','OK','OK');
uiputfile('*','delete box');
if size(file_list,1)==1
    file_list=gds.diagram;
elseif (index_selected==1)
    file_list=vertcat(file_list(index_selected+1:end,1));
elseif (index_selected==size(file_list,1))
    file_list=vertcat(file_list(1:index_selected-1,1));
else
    file_list=vertcat(file_list(1:index_selected-1,1),file_list(index_selected+1:end,1));
end
if isempty(strmatch('diagram',file_list))
    file_list{end+1} = gds.diagram; 
end
[status,msg]=mkdir(dir,'diagram');
gds.point = '';gds.type = '';gds.curve.new = '';
gds.open.figuur = 0; 
gds.open.numeric_fig = 0;
gds.open.continuer = 0;
gds.open.integrator = 0;
set(0,'ShowHiddenHandles','on');
set(MC.mainwindow.compute,'enable','off');
set(MC.mainwindow.window,'enable','off');
set(MC.mainwindow.Type,'enable','on');
set(MC.mainwindow.select_userfunctions,'enable','on');
%starter-window open?
close(MC.starter);
% continuer-window open?
close(MC.continuer); 
% numeric window open?
close(MC.numeric_fig); 
%2D-plot open
close(MC.D2); 
%3D-plot open
close(MC.D3);
close(MC.integrator); 
    close(MC.PRC); 
    gds.open.PRC = 0;
    %PRC-plot open
    close(MC.dPRC);
    gds.open.dPRC = 0;
load_matcont;
set(handles.listbox_diagram,'String',file_list);
cd(path_sys);cd ..,

%-------------------------------------------------------------------
function actions_rename_Callback(h,eventdata,handles,varargin)
global gds path_sys;
index_selected = get(handles.listbox_diagram,'Value');
file_list = get(handles.listbox_diagram,'String');
dir1 = file_list{index_selected};
prompt  = {'Enter the name of the diagram'};
title   = 'Rename';
lines= 1;
def     = {dir1};
answer  = inputdlg(prompt,title,lines,def);
if isempty(answer)
    return;
end
dir1=fullfile(path_sys,gds.system,dir1);
cd(dir1);d=dir;
cd(path_sys);cd ..;
str2 = char(answer);
dir2 = fullfile(path_sys,gds.system,str2);
[status,msg]=mkdir(path_sys,gds.system);
direc=fullfile(path_sys,gds.system);
[status,msg]=mkdir(direc,str2);
if (status==2)
    warndlg('This diagram already exists');
    return;
end
gds.diagram=str2;
for i=1:size(d,1)
    if ~d(i).isdir
        source=fullfile(dir1,d(i).name);
        dest=fullfile(dir2,d(i).name);
        copyfile(source,dest,'writable');
    end
end
file_list{end+1} = gds.diagram;
set(handles.listbox_diagram,'String',file_list);
cd(path_sys);cd ..,
% --------------------------------------------------------------------
function cancel_diagram_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.cancel_diagram.
global path_sys;
close(handles.diagramfig);
cd(path_sys);cd ..; 


% --------------------------------------------------------------------
function actions_new_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.cancel_diagram.
global path_sys gds MC FigPos;
file_list = get(handles.listbox_diagram,'String');
prompt  = {'Enter the name of the diagram'};
title   = 'Rename';
lines= 1;
def     = {''};
answer  = inputdlg(prompt,title,lines,def);
if isempty(answer)
    return;
end
str2 = char(answer);
dir = fullfile(path_sys,gds.system);
[status,msg]=mkdir(path_sys,gds.system);
[status,msg]=mkdir(dir,str2);
if (status==2)
    warndlg('This diagram already exists');
end
gds.diagram=str2;
file_list{end+1} = gds.diagram;
gds.point='';gds.type='';gds.curve.new='';
gds.open.figuur=0;
gds.open.numeric_fig=0;
gds.open.continuer=0;
gds.open.integrator=0;
load_matcont;
set(handles.listbox_diagram,'String',file_list);
cd(path_sys);cd ..,
set(MC.mainwindow.compute,'Enable','off');
set(MC.mainwindow.window,'enable','off');
set(MC.mainwindow.Type,'enable','on');
set(MC.mainwindow.select_userfunctions,'enable','on');
%starter-window open?
close(MC.starter);
% continuer-window open?
close(MC.continuer);
% numeric window open?
close(MC.numeric_fig);
%2D-plot open
close(MC.D2);
%3D-plot open
close(MC.D3);
%PRC-plot open
close(MC.PRC);
gds.open.PRC = 0;
%dPRC-plot open
close(MC.dPRC);
gds.open.dPRC = 0;
close(MC.integrator);