function varargout = editsystem(varargin)
% EDITSYSTEM Application M-file for editsystem.fig
%    FIG = EDITSYSTEM launch editsystem GUI.
%    EDITSYSTEM('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 05-Apr-2002 13:23:21
global path_sys oldgds gds driver_window
if ~isempty(gds.system)
    oldgds = gds;
end
if nargin <= 0  % LAUNCH GUI
    if nargin==0
         initial_dir = path_sys;
    elseif nargin==1 && exist(varargin{1},'dir')
        initial_dir = varargin{1};
    else 
        errordlg('Input argument must be a valid directory',...
            'Input Argument Error!')
        return
    end
    
    fig = openfig(mfilename,'reuse');

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
    %Populate the listbox
    load_listbox(initial_dir,handles)
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


%-----------------------------------------------------------------
function load_listbox(dir_path,handles)
dir_struct = dir(dir_path);
[sorted_names,sorted_index] = sortrows({dir_struct.name}');
handles.file_names = sorted_names;
handles.is_dir = [dir_struct.isdir];
handles.sorted_index = sorted_index;
guidata(handles.system_files,handles);
only_m_files(handles);


%-------------------------------------------------------------------
function only_m_files(handles)
global gds
number_files=length(handles.sorted_index);
j=1;val=j;
handles.file='';
for i=1:number_files
    file_list=handles.file_names;
    filename = file_list{i};
    if  ~handles.is_dir(handles.sorted_index(i))
        [path,name,ext] = fileparts(filename);
        if strcmp(ext,'.m')&& ~strcmp(filename,'standard.m')    
            handles.file{j}=filename;
            if strcmp(filename(1:end-2),gds.system)
                val=j;
            end
            j=j+1;
        end            
    end
end    
set(handles.files,'String',handles.file,'Value',val);

%-------------------------------------------------------------------
function del(handles)
global gds path_sys;
index_selected = get(handles.files,'Value');
file_list = get(handles.files,'String');
file  = file_list{index_selected};
file1 = fullfile(path_sys,file);
file2 = strcat(strrep(file1,'.m',''),'.mat');
dir1  = strrep(file1,'.m','');
if exist(file1,'file')&& exist(file2,'file')&& ~strcmp(strrep(file,'.m',''),gds.system)
    delete(file1); 
    delete(file2);
    str = sprintf('!rmdir %s /s/q',dir1);
    outp = evalc(str);
    str = sprintf('!rmdir -R %s',dir1);
    outp = evalc(str);
    rehash;
    h = guidata(handles.system_files);
    load_listbox(path_sys,h);
elseif size(file_list,1)>1
    i = 1;%file = file_list{1};
    while strcmp(strrep(file_list{i},'.m',''),gds.system)
        i = i+1;
    end 
    set(handles.files,'Value',i);
    edit_Callback([],[],handles,[]);
    delete(file1);  delete(file2);
    str = sprintf('!rmdir %s /s/q',dir1);
    outp = evalc(str);
    str = sprintf('!rmdir -R %s',dir1);
    outp = evalc(str);
    editsystem;
elseif size(file_list,1)==1
    close(handles.system_files)
    delete(file1);
    str = sprintf('!rmdir %s /s/q',dir1);
    outp = evalc(str);
    str = sprintf('!rmdir -R %s',dir1);
    outp = evalc(str);
    editsystem;       
    set(0,'ShowHiddenHandles','on');
    h = findobj('Type','figure','tag','system_files');
    set(0,'ShowHiddenhandles','off');
    edit_Callback(h,[],guidata(h),[]);
    delete(file2);
else
    errordlg('It isn''t possible to delete this system, maybe the system is currently loaded');
end
rehash;
%---------------------------------------------------------
function edit_Callback(h, eventdata, handles, varargin)
global path_sys gds hds MC FigPos;
% Stub for Callback of the uicontrol handles.files.
index_selected = get(handles.files,'Value');
file_list = get(handles.files,'String');
if isempty(file_list)%no systems
    close(handles.system_files)              
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
    %3D-plot open
    try
        close(MC.PRC);
    catch
        MC.PRC = [];
    end
    %3D-plot open
    try
        close(MC.dPRC);
    catch
        MC.dPRC = [];
    end
    close(MC.integrator);
    if strcmp(eventdata,'edit')
              systems('new');                 
    else
        systems('init');
        set(MC.mainwindow.compute,'enable','off');        %hal=findobj('Type','uimenu','Tag','window');
        set(MC.mainwindow.window,'enable','off');%        hal=findobj('Type','uimenu','Tag','Type');
        set(MC.mainwindow.Type,'enable','off');%        hal=findobj('Type','uimenu','Tag','select_userfunctions');
        set(MC.mainwindow.select_userfunctions,'enable','off');%        set(0,'ShowHiddenHandles','off');
        load_matcont;  
    end
else
    filename = file_list{index_selected};
    [path,name,ext] = fileparts(filename);    
    if strcmp(ext,'.m')
        try close(handles.system_files)              
             %starter-window open?
             if ~isempty(MC.starter), gds.open.figuur = 1;else gds.open.figuur = 0;end
             close(MC.starter);
             % continuer-window open?
             if ~isempty(MC.continuer), gds.open.continuer = 1;else gds.open.continuer = 0;end
             close(MC.continuer);
             % numeric window open?
             if ~isempty(MC.numeric_fig), gds.open.numeric_fig = 1;else gds.open.numeric_fig = 0;end
             close(MC.numeric_fig);
             %2D-plot open
             if ~isempty(MC.D2), gds.open.D2 = size(MC.D2);else gds.open.D2 = 0;end
             close(MC.D2);             
             %3D-plot open
             if ~isempty(MC.D3), gds.open.D3 = size(MC.D3);else gds.open.D3 = 0;end
             close(MC.D3);           
             %3D-plot open
             if ~isempty(MC.PRC), gds.open.PRC = size(MC.PRC);else gds.open.PRC = 0;end
         try
             close(MC.PRC); 
                gds.open.PRC = 0;
        catch
                MC.PRC = [];
        end
             %3D-plot open
             if ~isempty(MC.dPRC), gds.open.dPRC = size(MC.dPRC);else gds.open.dPRC = 0;end
         try
             close(MC.dPRC); 
                gds.open.dPRC = 0;
        catch
                MC.dPRC = [];
        end
             if ~isempty(MC.integrator),gds.open.integrator = 1;else gds.open.integrator = 0;end;
             close(MC.integrator);
             file = fullfile(path_sys,gds.system);
             save(file,'gds','FigPos');
             file = strrep(filename,'.m','');
             file = fullfile(path_sys,file);
             file = strcat(file,'.mat');
             load(file);
             if isfield(gds,'type') && ~isempty(gds.type)
                 feval(strcat('gui_',deblank(gds.type)));
             end
             if strcmp(eventdata,'edit')
                 systems;                 
             else
                 hds = [];
                 set(MC.mainwindow.compute,'enable','off');
                 set(MC.mainwindow.window,'enable','off');
                 set(MC.mainwindow.Type,'enable','on');
                 set(MC.mainwindow.select_userfunctions,'enable','on');
                 if gds.open.figuur==1;starter;end
                 if gds.open.continuer==1;continuer;end
                 if gds.open.numeric_fig==1;numeric;end
                 if gds.open.D2>0,for i=1:gds.open.D2,D2;end; end
                 if gds.open.D3>0,for i=1:gds.open.D3,plotD3;end; end
                 if gds.open.PRC>0,for i=1:gds.open.PRC,PRC_plot;end; end
                 if gds.open.dPRC>0,for i=1:gds.open.dPRC,dPRC_plot;end; end
                 if gds.open.integrator==1;integrator;end
                 load_matcont;
             end
         catch
             errordlg(lasterr,'File Type Error','modal');
             
         end
     end
 end



