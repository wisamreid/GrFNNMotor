function varargout = PRC_plot(varargin)
% PRC_window Application M-file
%    FIG = PRC_window launch PRC_window GUI.
%    PRC_window('callback_name', ...) invoke the named callback.

global gds MC FigPos;
if nargin == 0  % LAUNCH GUI
    fig = figure;
    nr=size(MC.PRC,2)+1;
%     set(fig,'KeyPressFcn','stop(''teststop'',gcbo)','Units','normalized','CloseRequestFcn','delete(gcf);'); %global MC; d = find(MC.PRC==gcf); MC.PRC(d) = []; 
    set(fig,'KeyPressFcn','stop(''teststop'',gcbo)','Units','normalized','CloseRequestFcn','PRC_plot(''closePRC'');'); %global MC; d = find(MC.PRC==gcf); MC.PRC(d) = []; 
    dat.nr=nr;dat.p2=[];dat.label=[];dat.traj=[];dat.labels=[];dat.trajs=[];
     set(fig,'color',[1 1 1],'UserData',dat);
    MC.PRC=[MC.PRC fig];
    nrfig=size(MC.PRC,2);
    for siz=1:nrfig
        set(MC.PRC(siz),'NumberTitle','off','Name',sprintf('Phase Response Curve:%d',siz));
    end
   	cla;
    % Generate a structure of handles to pass to callbacks, and store it.
    handles = guihandles(fig);
    guidata(fig, handles);
    if ~isfield(FigPos,'PRC') 
        FigPos.PRC= [];
    end
    if size(FigPos.PRC,2)< nr, FigPos.PRC{nr}=[0.5254 0.2826 0.4697 0.5677];end
    set(fig,'Position',FigPos.PRC{nr});
    set(gca,'XLIM',[0,1]);
    reset_callback(nrfig);
    if nargout > 0
		varargout{1} = fig;
    end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
		errordlg(lasterr,'error PRC');
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
if (isfield(gds.PRC(nr,1),'lim')&& ~isempty(gds.PRC(nr,1).lim))
    set(gca,'XLIM',[0,1]);
end
reset_callback(nr);

%----------------------------------------------------------------------
function redrawcurve(varargin)
global gds path_sys MC ; %lds;
if isnumeric(varargin{1}),nr=varargin{1};figure(MC.PRC(nr));
else d=get(gcf,'UserData');nr=d.nr;figure(gcf);end
file=fullfile(path_sys,gds.system);save(file,'gds');
file=strcat(gds.curve.new,'.mat');
file=fullfile(path_sys,gds.system,gds.diagram,file);
if exist(file,'file')
    load(file);
    string=feval(gds.gui_load_draw);
    if isempty(strmatch(char(gds.PRC(nr,1).type),string,'exact'))|isempty(strmatch(char(gds.PRC(nr,2).type),string,'exact'))
        PRC=[inf inf];
    else
        PRC = feval(gds.gui_make_ready,gds.PRC(nr,:),x);
    end
    feval(gds.gui_draw,2,PRC,file,'redraw');
else
    return;
end
%------------------------------------------------
function reset_callback(varargin)
global gds MC;
if (nargin==1)&&(isnumeric(varargin{1})),nr=varargin{1};h=get(MC.PRC(nr),'CurrentAxes');
else d = get(gcf,'Name');nr=str2double(strrep(d,'2Dplot:',''));h=get(gcf,'CurrentAxes');end
axes(h);d=axis;
p=get(h,'View');
xmode = get(h,'Xdir');
ymode = get(h,'YDir');
cla reset;
user=get(MC.PRC(nr),'UserData');
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
        PRC=feval(gds.gui_make_ready,gds.PRC(nr,:),x);
        feval(gds.gui_draw,2,PRC,filen,'redraw');
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
function closePRC(varargin)
global gds path_sys MC calculation_progress FigPos;%lds;
if calculation_progress~=0 
   return;
end;

d = find(MC.PRC==gcf);
user = get(gcf,'UserData');
if isfield(user,'nr')
    FigPos.PRC{user.nr} = get(gcf,'Position');
end
MC.PRC(d)=[];
if isempty(MC.PRC) || length(MC.PRC) == 0
    gds.open.PRC = 0;
    gds.options.PRC = 0;
    MC.PRC = [];
end
delete(gcf);
starter;