function varargout = dPRC_plot(varargin)
% dPRC_window Application M-file
%    FIG = dPRC_window launch dPRC_window GUI.
%    dPRC_window('callback_name', ...) invoke the named callback.

global gds MC FigPos;
if nargin == 0  % LAUNCH GUI
    fig = figure;
    nr=size(MC.dPRC,2)+1;
    set(fig,'KeyPressFcn','stop(''teststop'',gcbo)','Units','normalized','CloseRequestFcn','dPRC_plot(''closedPRC'');'); %global MC; d = find(MC.dPRC==gcf); MC.dPRC(d) = []; 
    dat.nr=nr;dat.p2=[];dat.label=[];dat.traj=[];dat.labels=[];dat.trajs=[];
     set(fig,'color',[1 1 1],'UserData',dat);
    MC.dPRC=[MC.dPRC fig];
    nrfig=size(MC.dPRC,2);
    for siz=1:nrfig
        set(MC.dPRC(siz),'NumberTitle','off','Name',sprintf('Derivative Phase Response Curve:%d',siz));
    end
   	cla;
    % Generate a structure of handles to pass to callbacks, and store it.
    handles = guihandles(fig);
    guidata(fig, handles);
    if ~isfield(FigPos,'dPRC') 
        FigPos.dPRC= [];
    end
    if size(FigPos.dPRC,2)< nr, FigPos.dPRC{nr}=[0.5254 0.2826 0.4697 0.5677];end
    set(fig,'Position',FigPos.dPRC{nr});
    set(gca,'XLIM',[0,1]);
    reset_callback(nrfig);
    if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
		errordlg(lasterr,'error dPRC');
        set(0,'ShowHiddenHandles','on');
        h=findobj('type','figure','Tag','stopcont');delete(h);
        set(0,'ShowHiddenHandles','off');   
	end

end
%--------------------------------------------------------------------
function OK(number)
global gds;
nr=number;
handles = guihandles(gcf);
guidata(gcf, handles);
if (isfield(gds.dPRC(nr,1),'lim')&& ~isempty(gds.dPRC(nr,1).lim))
    set(gca,'XLIM',[0,1]);
end
reset_callback(nr);

%----------------------------------------------------------------------
function redrawcurve(varargin)
global gds path_sys MC ; %lds;
if isnumeric(varargin{1}),nr=varargin{1};figure(MC.dPRC(nr));
else d=get(gcf,'UserData');nr=d.nr;figure(gcf);end
file=fullfile(path_sys,gds.system);save(file,'gds');
file=strcat(gds.curve.new,'.mat');
file=fullfile(path_sys,gds.system,gds.diagram,file);
if exist(file,'file')
    load(file);
    string=feval(gds.gui_load_draw);
    if isempty(strmatch(char(gds.dPRC(nr,1).type),string,'exact'))|isempty(strmatch(char(gds.dPRC(nr,2).type),string,'exact'))
        dPRC=[inf inf];
    else
        dPRC = feval(gds.gui_make_ready,gds.dPRC(nr,:),x);
    end
    feval(gds.gui_draw,2,dPRC,file,'redraw');
else
    return;
end
%------------------------------------------------
function reset_callback(varargin)
global gds MC;
if (nargin==1)&&(isnumeric(varargin{1})),nr=varargin{1};h=get(MC.dPRC(nr),'CurrentAxes');
else d = get(gcf,'Name');nr=str2double(strrep(d,'2Dplot:',''));h=get(gcf,'CurrentAxes');end
axes(h);d=axis;
p=get(h,'View');
xmode = get(h,'Xdir');
ymode = get(h,'YDir');
cla reset;
user=get(MC.dPRC(nr),'UserData');
set(h,'View',p,'XDir',xmode,'YDir',ymode);axis(d);

%----------------------------------------------------------------------
function redrawdiagram(varargin)
global gds path_sys ;%lds;
reset_callback;
d=get(gcf,'UserData');nr=d.nr;
oldgds = gds;
file=fullfile(path_sys,gds.system);save(file,'gds');
start_filen=gds.curve.new;
old_filen=gds.curve.old;
file=strcat(gds.curve.new,'.mat');
dir=fullfile(path_sys,gds.system,gds.diagram);
start_file=fullfile(dir,file);
cd(dir);w=what;
cd(path_sys);cd ..;
if isempty(w.mat)
    return;
else
    for i=1:size(w.mat,1)
        if strcmp(w.mat{i},'tempadwgyk_cycle.mat')
            continue;
        end
        filen=fullfile(path_sys,gds.system,gds.diagram,w.mat{i});
        load(filen);
        index = 1;
        if exist('cds','var')
            cds.options.UserfunctionsInfo=gds.options.UserfunctionsInfo;
            cds.options.Userfunctions=gds.options.Userfunctions;
            gds.options = cds.options;
        end
        gds.type = ctype; gds.point = point;
        feval(strcat('gui_',deblank(ctype)));
        gds.curve.new = strrep(w.mat{i},'.mat','');gds.curve.old = strrep(w.mat{i},'.mat',''); 
        file = fullfile(path_sys,gds.system);
        save(file,'gds');
        feval(gds.gui_load_point,index,x,'not_save',filen,1);      
        dPRC=feval(gds.gui_make_ready,gds.dPRC(nr,:),x);
        feval(gds.gui_draw,2,dPRC,filen,'redraw');
    end
end
if ~exist(start_file,'file')
    gds = oldgds;
     return;    
end
load(start_file); index = 1;
if exist('cds','var')
    cds.options.UserfunctionsInfo=gds.options.UserfunctionsInfo;
    cds.options.Userfunctions=gds.options.Userfunctions;
    gds.options = cds.options;
end
gds.type = ctype; gds.point = point;feval(strcat('gui_',deblank(ctype)));;
gds.curve.new = start_filen;gds.curve.old = old_filen; 
file = fullfile(path_sys,gds.system); save(file,'gds');
feval(gds.gui_load_point,index,x,'not_save',start_file,1);


%----------------------------------------------------------------------
function closedPRC(varargin)
global gds path_sys MC calculation_progress FigPos;%lds;
if calculation_progress~=0 
   return;
end;

d = find(MC.dPRC==gcf);
user = get(gcf,'UserData');
if isfield(user,'nr')
    FigPos.dPRC{user.nr} = get(gcf,'Position');
end
MC.dPRC(d)=[];
if isempty(MC.dPRC) || length(MC.dPRC) == 0
    gds.open.dPRC = 0;
    gds.options.dPRC = 0;
    MC.dPRC = [];
end
delete(gcf);
starter;
