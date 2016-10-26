function varargout = select_curve(varargin)
% EDITSYSTEM Application M-file for editsystem.fig
%    FIG = EDITSYSTEM launch editsystem GUI.
%    EDITSYSTEM('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 22-Oct-2001 12:50:00
global gds path_sys driver_window ;
if nargin <= 0  % LAUNCH GUI
    if nargin==0
        initial_dir=fullfile(path_sys,gds.system,gds.diagram);
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
        if ishandle(driver_window),delete(driver_window);end
        
	end
    
end

% --------------------------------------------------------------------
function select_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.files.
global gds path_sys MC;
initial_dir = fullfile(path_sys,gds.system,gds.diagram);
index_selected = get(handles.points,'Value');
file_list = get(handles.points,'String');
if isempty(file_list)||isempty(file_list{1})
    gds.curve.new='';gds.point='';gds.type='';
    set(MC.mainwindow.window,'enable','off');
    set(MC.mainwindow.compute,'enable','off');
    load_matcont;
    close(handles.select_curve);
    return
end
filename = file_list{index_selected};
filen = fullfile(initial_dir,filename);
load([filen,'.mat']); index = 1;
if exist('cds','var')
    cds.options.UserfunctionsInfo = gds.options.UserfunctionsInfo;
    cds.options.Userfunctions     = gds.options.Userfunctions;
    gds.options = contset;
    gds.options = cds.options;
elseif exist('option','var')
    gds.integrator.options = option;
end
gds.type = ctype;
if ~strcmp(ctype,'LC ')
    MC.PRC = [];
    gds.open.PRC = 0;
    MC.dPRC = [];
    gds.open.dPRC = 0;
end
gds.point = point;feval(strcat('gui_',deblank(ctype)));
gds.curve.new = filename;gds.curve.old = filename; 
file = fullfile(path_sys,gds.system);
save(file,'gds');
feval(gds.gui_load_point,index,x,'save',filen,1);
list = get(MC.mainwindow.initial_point,'children');
for i = 1:length(list)
    tag = get(list(i),'Tag'); [tag1,tag2] = strtok(tag,'_');
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
load_matcont;
if exist('cds','var');continuer;else integrator;end
starter;%load the appropriate continuer and starter window (see starter)
set(MC.mainwindow.extend,'enable','on');
set(MC.mainwindow.forward,'enable','on');
set(MC.mainwindow.backward,'enable','on');
if ~isempty(MC.numeric_fig), numeric;end
close(handles.select_curve);

%-----------------------------------------------------------------
function load_listbox(dir_path,handles)
global path_sys gds;
val = 1; string = '';
if (exist(dir_path,'dir')~=7)
    string = '';
else
    cd (dir_path);
    w = what;sorted_names = sort(w.mat);
    j = strmatch('tempadwgyk_cycle.mat',sorted_names,'exact');
    if ~isempty(j),sorted_names(j)=[];end
    number_files = length(sorted_names);
    if number_files>0
        for i = 1:number_files
            sn = sorted_names{i};
            if size(sn,1) > 1
                sn = sn(1,:);
            end
            string{i,1} = strrep(sn,'.mat','');
            if strcmp(string{i,1},gds.curve.new)
                val = i;
            end
        end
    end
end
guidata(handles.select_curve,handles);
set(handles.points,'String',string,'Value',val);
cd(path_sys);cd ..;
%-------------------------------------------------------------------
function cancel_Callback(h, eventdata, handles, varargin)
global path_sys gds
file_list = get(handles.points,'String');
filemat = strcat(gds.curve.new,'.mat');
filemat = fullfile(path_sys,gds.system,gds.diagram,filemat);
if size(file_list,1)==1 && isempty(file_list{1})
    file_list=[];
end
if ~exist(filemat,'file')&& size(file_list,1)>=1
    i=1;
    while strcmp(file_list{i},gds.curve.new)
        i=i+1;
    end 
    set(handles.points,'Value',i);
    select_Callback([],[],handles,[]);
elseif ~exist(filemat,'file') && isempty(file_list)
    h=guidata(handles.select_curve);
    select_Callback(handles.select_curve,[],h,[]);
else
    close(handles.select_curve);
end
cd(path_sys);cd ..; 

%-------------------------------------------------------------------
function delete_callback(h,eventdata,handles,varargin)
global gds path_sys;
index_selected = get(handles.points,'Value');
file_list = get(handles.points,'String');
file    = file_list{index_selected};
filemat = strcat(file,'.mat');
filemat = fullfile(path_sys,gds.system,gds.diagram,filemat);
initial_dir = fullfile(path_sys,gds.system,gds.diagram);
if exist(filemat,'file')&& ~strcmp(file,gds.curve.new)   
    delete(filemat);
    rehash;
    h=guidata(handles.select_curve);
    load_listbox(initial_dir,h);
elseif size(file_list,1)>1
    i=1;
    while strcmp(file_list{i},gds.curve.new)
        i=i+1;
    end 
    set(handles.points,'Value',i);
    delete(filemat);
    rehash;
    h=guidata(handles.select_curve);
    load_listbox(initial_dir,h);
elseif size(file_list,1)==1
    delete(filemat);
    rehash;
    h=guidata(handles.select_curve);
    load_listbox(initial_dir,h);
else
    warndlg('You have to select a curve')
    return
end

%-------------------------------------------------------------------
function rename_callback(h,eventdata,handles,varargin)
global gds path_sys;
index_selected = get(handles.points,'Value');
file_list = get(handles.points,'String');
file = file_list{index_selected};
num = str2num(strtok(file,')'));
if isempty(num) 
    prompt  = {'Enter the name of the curve'};
    title   = 'Rename';
    lines= 1;
    def     = {file};
    answer  = inputdlg(prompt,title,lines,def);
    if isempty(answer)
        return;
    end
    str2 = char(strcat(answer,'.mat'));
    dir2 = fullfile(path_sys,gds.system,gds.diagram,str2);
    str1 = strcat(file,'.mat');
    file_list{index_selected} = strrep(str2,'.mat','');
    dir1 = fullfile(path_sys,gds.system,gds.diagram,str1);
    copyfile(dir1,dir2,'writable');
    delete(dir1);
    set(handles.points,'String',file_list);
else
    warndlg('You have to select a curve')
    return
end
