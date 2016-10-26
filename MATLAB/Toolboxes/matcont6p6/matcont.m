function varargout = matcont(varargin)
% MATCONT Application M-file for matcont.fig
%    FIG = MATCONT launch matcont GUI.
%    MATCONT('callback_name', ...) invoke the named callback.


global gds path_sys sys calculation_progress driver_window MC FigPos;

%add the path to the continuer
warning off
string = get(0,'defaultuicontrolfontname');
% set(0,'defaultuicontrolfontsize',12);
set(0,'fixedwidthfontname',string); 
if nargin == 0  % LAUNCH GUI
    [list,val] = spparms;
    spparms('default');
    for i=1:length(list(:,1))
        if strcmp(list(i,:),'umfpack')
%             spparms('umfpack',0); %switch umfpack off, use v4solver
        end
    end

    addpath(cd);
    addpath([cd '/Continuer/']);
    addpath([cd '/Equilibrium/']);
    addpath([cd '/LimitCycle/']);
    addpath([cd '/PeriodDoubling/']);
    addpath([cd '/Systems/']);
    addpath([cd '/GUI/']);
    addpath([cd '/LimitPoint/']);
    addpath([cd '/Hopf/']);
    addpath([cd '/LimitPointCycle/']);
    addpath([cd '/NeimarkSacker/']);
    addpath([cd '/BranchPoint/']);
    addpath([cd '/BranchPointCycle/']);
    addpath([cd '/Homoclinic/']);
    addpath([cd '/HomoclinicSaddleNode/']);
    addpath([cd '/ConnectionHom/']);
    addpath([cd '/ConnectionHSN/']);
    addpath([cd '/HomotopySaddle/']);
    addpath([cd '/HomotopySaddleNode/']);
    addpath([cd '/HomotopyHet/']);
    addpath([cd '/ConnectionHet/']);
    addpath([cd '/Heteroclinic/']);
    addpath([cd '/MultilinearForms/']);
    addpath([cd '/Help/']);
    addpath([cd '/LimitCycleCodim2/']);
    addpath([cd '/SBML/']);
    
%    addpath([cd '/Help/']);
    % Find and go to correct directory (class directory)
    p = mfilename('fullpath');
    p = p(1:length(p)-length(mfilename));
    p = strcat(p,'/LimitCycle');
    curdir = cd;
    cd(p);
    
    % Compile the c-files (optimized)
   if ~(exist(strcat('BVP_LC_jac.',mexext),'file'))
   if ~isempty(regexp(mexext,'64','match'))
        mex -largeArrayDims -O BVP_LC_jac.c;
        mex -largeArrayDims -O BVP_PD_jac.c;
        mex -largeArrayDims -O BVP_BPC_jacC.c;
        mex -largeArrayDims -O BVP_BPC_jacCC.c;
        mex -largeArrayDims -O BVP_LPC_jac.c;
        mex -largeArrayDims -O BVP_NS_jac.c;
        mex -largeArrayDims -O BVP_LCX_jac.c;
   else
        mex -O BVP_LC_jac.c;
        mex -O BVP_PD_jac.c;
        mex -O BVP_BPC_jacC.c;
        mex -O BVP_BPC_jacCC.c;
        mex -O BVP_LPC_jac.c;
        mex -O BVP_NS_jac.c;
        mex -O BVP_LCX_jac.c;
  end
end

    
    % Return to directory we started in
    cd (curdir);
    
    
    MC=[];driver_window=[];calculation_progress=0;
    MC.mainwindow=[];MC.starter=[];MC.continuer=[];
    MC.integrator=[];MC.numeric_fig=[];MC.D2=[];MC.PRC=[];MC.dPRC=[];
    MC.D3=[];
    fig = openfig(mfilename,'reuse','invisible');
    handles = guihandles(fig); 
    guidata(fig, handles);
    MC.mainwindow = handles;
    %set the handle of the window visible so you
    %can change the menu of it
	set(fig,'Color',[1 1 1]);
    path_sys = fullfile(pwd,'Systems');
    gds = [];gds.system = '';x=[];v=[];s=[];eds=[];
    %is Matcont already been used?    
    session = fullfile(path_sys,'session.mat');
    %case yes: load session file and extract file from it
    if (exist(session,'file')==2) 
        load(session);       
        file = fullfile(path_sys,sys.file);
    end
    % if the previous calls worked, then test if the file exists
    if exist('file','var')&&exist(file,'file')
        load(file);
        if isfield(gds,'type') & ~isempty(gds.type)
            feval(strcat('gui_',deblank(gds.type)))
        end
        if isempty(FigPos) 
            FigPos=[];FigPos.mainwindow=[5 5 10 10];FigPos.starter=[90.0 3 50 18];FigPos.continuer=[150 2.995 48.167 31.438];
            FigPos.integrator=[150 2.995 43 26.692];FigPos.numeric_fig=[120 3. 47.0000  15];FigPos.D2{1}=[0.5254 0.2826 0.4697 0.5677];
            FigPos.D3{1}=[0.004 0.335 0.459 0.439];FigPos.PRC{1}=[0.004 0.335 0.459 0.439];FigPos.dPRC{1}=[0.004 0.335 0.459 0.439];
        end
        if ~isempty(FigPos.mainwindow)
            set(fig,'Position',FigPos.mainwindow);
        end
        movegui(fig,'onscreen');
        set(fig,'Name','MatCont');
        set(fig, 'Visible','on');
        if gds.open.figuur==1;starter;end
        if gds.open.continuer==1;continuer;end
        if gds.open.numeric_fig==1;numeric;end
        if gds.open.D2>0,for i=1:gds.open.D2,D2;end; end
        if gds.open.D3>0,for i=1:gds.open.D3,plotD3;end; end
        if gds.open.PRC>0,for i=1:gds.open.PRC,PRC_plot;end; end
        if gds.open.dPRC>0,for i=1:gds.open.dPRC,dPRC_plot;end; end
        if gds.open.integrator==1;integrator;end
        if isempty(gds.type)
            set(MC.mainwindow.window,'Enable','off');
        end
    else %first matcont runs
        set(MC.mainwindow.compute,'enable','off');
        set(MC.mainwindow.window,'enable','off');
        set(MC.mainwindow.Type,'enable','off');  
        set(MC.mainwindow.select_userfunctions,'enable','off');
        sys.gui.pausenever=0;
        sys.gui.pauseeachpoint=0;
        sys.gui.pausespecial=1;
        sys.gui.ncurves=2;
        sys.gui.plot_points=1;
        if isempty(FigPos) 
            FigPos=[];FigPos.mainwindow=[];FigPos.starter=[90.0 3 50 18];FigPos.continuer=[150 0.563 48.167 31.438];
            FigPos.integrator=[150 2.995 43 26.692];FigPos.numeric_fig=[120 3. 47.0000  15];FigPos.D2{1}=[0.5254 0.2826 0.4697 0.5677];
            FigPos.D3{1}=[0.004 0.335 0.459 0.439];FigPos.PRC{1}=[0.004 0.335 0.459 0.439];FigPos.dPRC{1}=[0.004 0.335 0.459 0.439];
        end
        movegui(fig,'onscreen');
        set(fig,'Name','MatCont');
        set(fig, 'Visible','on');
    end
    if ~isfield(sys.gui,'plot_points'),sys.gui.plot_points=1;end
    if ~isfield(sys.gui,'ncurves'),sys.gui.ncurves=2;end
    set(MC.mainwindow.options_window,'Userdata',sys.gui.plot_points);
    % Generate a structure of handles to pass to callbacks, and store it. 
    if ~isempty(gds.system)
        load_matcont;
    end
    if nargout > 0
		varargout{1} = fig;
    end
elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

        if isempty(FigPos) 
            FigPos=[];FigPos.mainwindow=[];FigPos.starter=[90.0 3 50 18];FigPos.continuer=[150 0.563 48.167 31.438];
            FigPos.integrator=[150 2.995 43 26.692];FigPos.numeric_fig=[120 3. 47.0000  15];FigPos.D2{1}=[0.5254 0.2826 0.4697 0.5677];
            FigPos.D3{1}=[0.004 0.335 0.459 0.439];FigPos.PRC{1}=[0.004 0.335 0.459 0.439];FigPos.dPRC{1}=[0.004 0.335 0.459 0.439];
        end

              	try
        varargin{1} = str2func(varargin{1});
    	[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
              	catch
		errordlg(lasterr,'error matcont');
        calculation_progress=0;
        if ishandle(driver_window),delete(driver_window);end
        if ishandle(MC.mainwindow.window), set(MC.mainwindow.window,'Enable','on');end
        if ishandle(MC.mainwindow.select),set(MC.mainwindow.select,'enable','on');end
        if ishandle(MC.mainwindow.options),set(MC.mainwindow.options,'enable','on'); end
        if ishandle(MC.mainwindow.compute),set(MC.mainwindow.compute,'enable','on');end
        if ishandle(MC.mainwindow.Type),set(MC.mainwindow.Type,'enable','on');end
        if ishandle(MC.mainwindow.extend),set(MC.mainwindow.extend,'enable','on');end
        if ishandle(MC.mainwindow.forward),set(MC.mainwindow.forward,'enable','on');end
        if ishandle(MC.mainwindow.backward),set(MC.mainwindow.backward,'enable','on'); end    
              	end

end

%----------------------------------------------------------------
function exit_callback(h,eventdata,handles,varargin)
%save the last used system and current settings
global gds sys path_sys MC FigPos;
selection = questdlg('Exit  MatCont?',...
                     'Close MatCont',...
                     'Yes','No','Yes');
switch selection
case 'No'
     return
end
shh=get(0,'ShowHiddenHandles');
set(0,'ShowHiddenHandles','on');
currFig=get(0,'CurrentFigure');
%starter-window open?
if ~isempty(MC.starter)
    close(MC.starter);
    gds.open.figuur=1;
else gds.open.figuur=0;
end
% continuer-window open?
if ~isempty(MC.continuer)
    close(MC.continuer);
    gds.open.continuer=1;
else gds.open.continuer=0;
end
% numeric window open?
if ~isempty(MC.numeric_fig)
    close(MC.numeric_fig);
    gds.open.numeric_fig=1;
else gds.open.numeric_fig=0;
end
%2D-plot open
if size(MC.D2,2)>0
    gds.open.D2 = size(MC.D2,2);
    for k=1:gds.open.D2,        
        close(MC.D2(1));
    end
else gds.open.D2=0;
end
%3D-plot open
if size(MC.D3,2)>0
    gds.open.D3=size(MC.D3,2);
    for k=1:gds.open.D3,
        close(MC.D3(1));
    end
else gds.open.D3=0;
end
%PRC-plot open
if size(MC.PRC,2)>0
    gds.open.PRC = size(MC.PRC,2);
    for k=1:gds.open.PRC       
        figure(MC.PRC(1));
        close(MC.PRC(1));
    end
else gds.open.PRC=0;
end
%dPRC-plot open
if size(MC.dPRC,2)>0
    gds.open.dPRC = size(MC.dPRC,2);
    for k=1:gds.open.dPRC,       
        figure(MC.dPRC(1));     
        close(MC.dPRC(1));
    end
else gds.open.dPRC=0;
end
%h = findobj('type','figure','tag','integrator');delete(h);
if ~isempty(MC.integrator)
    close(MC.integrator);
    gds.open.integrator=1;
else gds.open.integrator=0;
end
FigPos.mainwindow = get(currFig,'Position');
set(0,'ShowHiddenHandles',shh); 
if ~isfield(gds,'system') || isempty(gds.system)
    delete(currFig);
    close all hidden;
    return    
else
    file = fullfile(path_sys,gds.system);
    save(file,'gds','FigPos');
    sys.file = strcat(gds.system,'.mat');
   % close all hidden;
    session = fullfile(path_sys,'session.mat');
    save(session,'sys');
end
clear FigPos
delete(currFig);

%----------------------------------------------------------------
function forward_callback(h,eventdata,handles,varargin)
%if you press the menu 'compute->forward'at the matcont-window
%the continuer will start and the direction will be set to forward
global gds cds path_sys MC calculation_progress;
gds.options.Backward=0;
make_curvename;
load_matcont;
set(MC.mainwindow.window,'enable','off');
set(MC.mainwindow.select,'enable','off');
set(MC.mainwindow.options,'enable','off');
set(MC.mainwindow.compute,'enable','off');
set(MC.mainwindow.Type,'enable','off');
if ~isempty(MC.numeric_fig),numeric;drawnow;end
calculation_progress = 1;
start_cont;
calculation_progress = 0;

%-----------------------------------------------------------------
function backward_callback(h,eventdata,handles,varargin)
%if you press the menu 'compute->backward'at the matcont-window
%the continuer will start and the direction will be set to backward
global cds gds path_sys MC calculation_progress;
gds.options.Backward = 1;
make_curvename;
load_matcont;
set(MC.mainwindow.window,'enable','off');
set(MC.mainwindow.select,'enable','off');
set(MC.mainwindow.options,'enable','off');
set(MC.mainwindow.compute,'enable','off');
set(MC.mainwindow.Type,'enable','off');
if ~isempty(MC.numeric_fig),numeric;drawnow;end
calculation_progress = 1;
start_cont;
calculation_progress = 0;

%-----------------------------------------------------------------
function extend_callback(h,eventdata,handles,varargin)
%if you press the menu 'compute->extend'at the matcont-window
%the continuer will start and the direction will be the same as before
global MC calculation_progress;
load_matcont;
set(MC.mainwindow.window,'enable','off');
set(MC.mainwindow.select,'enable','off');
set(MC.mainwindow.options,'enable','off');
set(MC.mainwindow.compute,'enable','off');
set(MC.mainwindow.Type,'enable','off');
calculation_progress = 1;
start_cont('extend');
calculation_progress = 0;

%------------------------------------------------------------------
function point_callback(handles)
global gds MC;
tag = get(gcbo,'Tag');
[point,str] = strtok(tag,'_');
if ~isempty(str)
    if ~isempty(gds.type)
        type = strcat('_',deblank(gds.type));
    else
        type='';
    end
    if isempty(findstr(str,type))
        type = strtok(str,'_');
    else
        type = strtok(type,'_');
    end
    if ~strcmp(type,'LC')
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
    feval(strcat('gui_',type));
    feval(gds.gui_point,tag);
else
    str = sprintf('It isn''t provided to start a continuation from a %s point',point);
    errordlg(str);
end

%------------------------------------------------------------------
function type_callback(varargin)
global gds MC
newtag = get(gcbo,'Tag');
if ~strcmp(newtag,'LC ')
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
if strcmp(newtag,'NS2')
    gds.whichNS = 2;
    newtag = 'NS ';
elseif strcmp(newtag,'NS')
    gds.whichNS = 1;
end
type =strcat('gui_',newtag);
feval(type);
feval(gds.gui_type);


%------------------------------------------------------------------
function make_curvename
%makes a non-existing name for the new curve
global gds path_sys sys MC;
str = strcat(gds.point,'_',gds.type,'(1).mat');
[status,msg] = mkdir(path_sys,gds.system);
file=fullfile(path_sys,gds.system,gds.diagram,str);
i=1;        
while (exist(file,'file')~=0)&&(i<sys.gui.ncurves)
      d(i) = dir(file);
      i = i+1;
      str = strrep(str,strcat('(',num2str(i-1)),strcat('(',num2str(i)));
      file = fullfile(path_sys,gds.system,gds.diagram,str);
end

if (i==sys.gui.ncurves)&&(exist(file,'file')~=0)
    d(i) = dir(file);
    % get the datenumber from the file, first try the dir.datenum field,
    % this works in more locales and is the correct way to get it, but
    % versions prior to MATLAB R2007a did not contain that field, so in
    % that case, try an error prone method of getting the datenum
    try
        da = vertcat(d.datenum);
    catch
        da=datenum(vertcat(d(:).date), 'dd-mmm-yyyy HH:MM:SS');
    end
    j = find(da==min(da));
    j = j(1);
    str=strrep(str,num2str(i),num2str(j));
    file=fullfile(path_sys,gds.system,gds.diagram,str);
end
gds.curve.new=strrep(str,'.mat','');
file=fullfile(path_sys,gds.system);
save(file,'gds');
set(MC.mainwindow.extend,'enable','off');

%------------------------------------------------------------------
function make_curve
%makes a non-existing name for the new curve
global gds path_sys sys MC;
str = strcat(gds.point,'_',gds.type,'(1).mat');
[status,msg] = mkdir(path_sys,gds.system);
file=fullfile(path_sys,gds.system,gds.diagram,str);
i=1;        
while (exist(file,'file')~=0)&&(i<sys.gui.ncurves)
      d(i) = dir(file);
      i = i+1;
      str = strrep(str,strcat('(',num2str(i-1)),strcat('(',num2str(i)));
      file = fullfile(path_sys,gds.system,gds.diagram,str);
end
if (i==sys.gui.ncurves)&&(exist(file,'file')~=0)
    d(i) = dir(file);
    % get the datenumber from the file, first try the dir.datenum field,
    % this works in more locales and is the correct way to get it, but
    % versions prior to MATLAB R2007a did not contain that field, so in
    % that case, try an error prone method of getting the datenum
    try
        da = vertcat(d.datenum);
    catch
        da=datenum(vertcat(d(:).date), 'dd-mmm-yyyy HH:MM:SS');
    end
    j = find(da==min(da));
    j = j(1);
    str=strrep(str,num2str(i),num2str(j));
    file=fullfile(path_sys,gds.system,gds.diagram,str);
    delete(file);
    struct1 = contget(gds.options,'UserfunctionsInfo',[]);
    struct2 = contget(gds.options,'Userfunctions',0);
    gds.options = contset;
    gds.options = contset(gds.options,'Singularities',1);
    gds.options = contset(gds.options,'UserfunctionsInfo',struct1);
    gds.options = contset(gds.options,'Userfunctions',struct2);
    gds.options.IgnoreSingularity=[];
    gds.options.Multipliers = 1; 
    gds.options.Eigenvalues = 1;
    gds.options.Locators = [];
end
gds.curve.new=strrep(str,'.mat','');
file=fullfile(path_sys,gds.system);
save(file,'gds');
set(MC.mainwindow.extend,'enable','off');

%------------------------------------------------------------------
function start_cont(varargin)
global gds path_sys sys MC driver_window calculation_progress;
file=fullfile(path_sys,gds.system);
save(file,'gds');
calculation_progress=1;
if ~isempty(gds.type)
    if isempty(varargin)
        feval(gds.gui_start_cont);
    else
       feval(gds.gui_start_cont,'extend');  
   end
else
    [x,v,s]=cont(gds.system,x0,[],gds.options);
end
set(MC.mainwindow.window,'enable','on');
set(MC.mainwindow.select,'enable','on');
set(MC.mainwindow.options,'enable','on');
set(MC.mainwindow.compute,'enable','on');
set(MC.mainwindow.extend,'enable','on');
set(MC.mainwindow.forward,'enable','on');
set(MC.mainwindow.backward,'enable','on');
set(MC.mainwindow.Type,'enable','on');


% --------------------------------------------------------------------
function select_exit_Callback(hObject, eventdata, handles)
% hObject    handle to select_exit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


