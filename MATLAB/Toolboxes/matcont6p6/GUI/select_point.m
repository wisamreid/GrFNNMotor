function varargout = select_point(varargin)
% EDITSYSTEM Application M-file for editsystem.fig
%    FIG = EDITSYSTEM launch editsystem GUI.
%    EDITSYSTEM('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 22-Oct-2001 12:50:00
global gds path_sys driver_window 
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
        if ishandle(driver_window),delete(driver_window); end
    end
    
end

%---------------------------------------------------------------------
function points_Callback(h,eventdata,handles,varargin)
global path_sys gds MC;
if strcmp(get(handles.select_point,'SelectionType'),'open');
    if (isempty(MC.D2) && isempty(MC.D3))
        str = sprintf('1)load the selected curve\n2)open window and set the attributes');
        errordlg(str,'show point on curve');
        return;
    else
        initial_dir = fullfile(path_sys,gds.system,gds.diagram);
        cd(initial_dir);
        index_selected = get(handles.points,'Value');
        file_list = get(handles.points,'String');
        points = file_list{index_selected};
        num = str2double(strtok(points,')'));
        if isempty(num)
            return
        end
        index_selected = index_selected-num;
        filename = file_list{index_selected};
        if (~strcmp(filename,gds.curve.new))
            str = sprintf('1)load the selected curve\n2)open window and set the attributes');
            errordlg(str,'show point on curve');   
            return;
        end
        h = []; f = [];
        load([filename,'.mat']);s2=s;
        num   = str2double(strtok(points,')'));
        index = s2(num).index;
        feval(strcat('gui_',deblank(gds.type)));
        if ~strcmp(deblank(gds.type),'LC')
            if ~isempty(MC.PRC)
                close(MC.PRC)
            end
            MC.PRC = [];
            gds.open.PRC = 0;
            if ~isempty(MC.dPRC)
                close(MC.dPRC)
            end
            MC.dPRC = [];
            gds.open.dPRC = 0;
        end
        feval(gds.gui_load_point,index,x,'not_save',filename,num);
        if  size(MC.D2,2)>0
            for siz = 1:size(MC.D2,2)
                dat   = get(MC.D2(siz),'UserData');nr=dat.nr;
                plot2 = feval(gds.gui_make_ready,gds.plot2(nr,:),x);
                figure(MC.D2(siz));
                s1 = s2(num);
                d = axis; D2('reset_callback');
                D2('redrawcurve',guidata(gcbo));
                feval(gds.gui_draw,2,plot2,x,'select',s1,h,f);
            end
        end
        load([filename,'.mat'])
        s2=s;
        num   = str2double(strtok(points,')'));
        index = s2(num).index;
        if size(MC.D3,2)>0
            for siz=1:size(MC.D3,2)
                dat=get(MC.D3(siz),'UserData');nr=dat.nr;
                plo=feval(gds.gui_make_ready,gds.plot3(nr,:),x);
                figure(MC.D3(siz));s1 = s2(num); d=axis;plotD3('reset_callback');
                plotD3('redrawcurve',guidata(gcbo));
                feval(gds.gui_draw,3,plo,x,'select',s1,h,f);          
            end
        end
        load([filename,'.mat']);s2=s;
        num   = str2double(strtok(points,')'));
        index = s2(num).index;
        if size(MC.numeric_fig,2)>0
            if strcmp(ctype,'DO ') || strcmp(ctype,'O ')
                d = get(MC.numeric_fig,'Userdata');
                d.label=feval(gds.gui_numeric_label,[]);
                set(MC.numeric_fig,'Userdata',d);
                figure(MC.numeric_fig);s1 = s2(num);% d=axis;plotD('reset_callback');
                feval(gds.gui_draw,4,t,x,'select',s1);          
            else
                numsing=sing_numeric(cds.options.IgnoreSingularity);
                d = get(MC.numeric_fig,'Userdata');
                d.label=feval(gds.gui_numeric_label,numsing);
                set(MC.numeric_fig,'Userdata',d);
                figure(MC.numeric_fig);s1 = s2(num);% d=axis;plotD('reset_callback');
                feval(gds.gui_draw,4,[],x,'select',s1,h,f);          
            end
        end
        feval(gds.gui_load_point,1,x,'not_save',filename,1);
    end
end
cd(path_sys);cd ..;

% --------------------------------------------------------------------
function select_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.files.
global gds path_sys MC;
initial_dir = fullfile(path_sys,gds.system,gds.diagram);
index_selected = get(handles.points,'Value');
file_list = get(handles.points,'String');
if isempty(file_list), close(handles.select_point);return; end
points = file_list{index_selected};
num = str2double(strtok(points,')'));
if isempty(num)
    return
end
index_selected = index_selected-num;
filename = file_list{index_selected};
filename = strcat(filename,'.mat');
file = fullfile(initial_dir,filename);
load(file);feval(strcat('gui_',deblank(ctype)));
        if ~strcmp(deblank(ctype),'LC')
            if ~isempty(MC.PRC)
                close(MC.PRC)
            end
            MC.PRC = [];
            gds.open.PRC = 0;
            if ~isempty(MC.dPRC)
                close(MC.dPRC)
            end
            MC.dPRC = [];
            gds.open.dPRC = 0;
        end
num = str2double(strtok(points,')'));
index = s(num).index;
            if size(filename,1) > 1
                sn = sn(filename,:);
            end
gds.curve.old = strrep(filename,'.mat','');
feval(gds.gui_load_point,index,x,'save',file,num);
if strcmp(s(num).label,'00')||strcmp(s(num).label,'99')
    s(num).label = ctype;
end
file = fullfile(path_sys,gds.system);
save(file,'gds');
list=get(MC.mainwindow.initial_point,'children');
if isfield(MC.mainwindow,'codim2')
    list2=get(MC.mainwindow.codim2,'children');
else
    list2 = [];
end
list=[list;list2];
tag = get(list,'Tag');
slabel = deblank(s(num).label);
if strmatch('BP',s(num).label),slabel='BP';end
if strmatch('BPC',s(num).label),slabel='BPC';end
for i = 1:length(list)
    [tag1,tag2] = strtok(tag{i},'_');
    if strcmp(slabel,tag1)
        type = strtok(tag2,'_');
        tag{i}=strcat(s(num).label,'_',tag2);
        break
    else
        type = deblank(ctype);
    end
end
if isempty(type)
    feval(strcat('gui_',deblank(ctype)));
    feval(gds.gui_point,tag{i});
else
    feval(strcat('gui_',type));
    feval(gds.gui_point,tag{i});
end
load_matcont;
if size(MC.numeric_fig)==1; numeric;end
close(handles.select_point);

%-----------------------------------------------------------------
function load_listbox(dir_path,handles)
global path_sys gds;
if exist(dir_path,'dir')
    cd (dir_path);
    w = what;
else
    string ='';val =1;
end
%if (exist(gds.system,'dir')~=7)|size(w.mat,1)==0
if size(w.mat,1)==0
    string = ''; val  = 1;
else
    sorted_names = sort(w.mat);
    j = strmatch('tempadwgyk_cycle.mat',sorted_names,'exact');
    if ~isempty(j),sorted_names(j)=[];end
    number_files = length(sorted_names);
    if number_files>0
        p = 1; val = p;
        for i=1:number_files
            sn = sorted_names{i};
            if size(sn,1) > 1
                sn = sn(1,:);
            end
            string{p,1}=strrep(sn,'.mat','');
            if strcmp(string{p,1},gds.curve.new)
                val = p;
            end
            load(sorted_names{i});
            p = p+1;
            j = length(s);
            for k=1:j
                string{p,1} = sprintf('%s) %s:%c%s',num2str(k),s(k).label,char(32),s(k).msg);
                p = p+1;
            end
        end
    end
end
guidata(handles.select_point,handles);
set(handles.points,'String',string,'Value',val);
cd(path_sys);cd ..;
%-------------------------------------------------------------------
function cancel_Callback(h, eventdata, handles, varargin)
global path_sys 
close(handles.select_point);
cd(path_sys);cd ..; 

%-------------------------------------------------------------------
function delete_callback(h,eventdata,handles,varargin)
global gds path_sys;
index_selected = get(handles.points,'Value');
file_list = get(handles.points,'String');
file = file_list{index_selected};
file = strcat(file,'.mat');
file = fullfile(path_sys,gds.system,gds.diagram,file);
initial_dir = fullfile(path_sys,gds.system,gds.diagram);
if exist(file,'file')&& ~strcmp(file,gds.curve.new)   
    delete(file);
    rehash;
    h=guidata(handles.select_point);
    load_listbox(initial_dir,h);
elseif size(file_list,1)>1
    i=1;
    while strcmp(file_list{i},gds.curve.new)
        i=i+1;
    end 
    set(handles.points,'Value',i);
    delete(file);
    rehash;
    h=guidata(handles.select_point);
    load_listbox(initial_dir,h);
elseif size(file_list,1)==1
    delete(file);
    rehash;
    h=guidata(handles.select_point);
    load_listbox(initial_dir,h);
else
    warndlg('You have to select a curve')
    return
end

%-------------------------------------------------------------------
function rename_callback(h,eventdata,handles,varargin)
global gds path_sys;
index_selected=get(handles.points,'Value');
file_list=get(handles.points,'String');
file=file_list{index_selected};
num=str2double(strtok(file,')'));
if isempty(num) 
    prompt  = {'Enter the name of the curve'};
    title   = 'Rename';
    lines= 1;
    def     = {file};
    answer  = inputdlg(prompt,title,lines,def);
    if isempty(answer),return;end
    str2=char(strcat(answer,'.mat'));
    dir2=fullfile(path_sys,gds.system,gds.diagram,str2);
    str1=strcat(file,'.mat');
    file_list{index_selected}=strrep(str2,'.mat','');
    dir1=fullfile(path_sys,gds.system,gds.diagram,str1);
    copyfile(dir1,dir2,'writable');
    delete(dir1);
    set(handles.points,'String',file_list);
else
    warndlg('You have to select a curve')
    return
end

