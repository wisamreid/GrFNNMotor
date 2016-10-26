function varargout = D2(varargin)
% D2 Application M-file for D2.fig
%    FIG = D2 launch D2 GUI.
%    D2('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 25-Sep-2002 13:34:14
global gds MC FigPos;
if nargin == 0  % LAUNCH GUI
    fig = openfig(mfilename,'new','invisible');
    set(0,'ShowHiddenHandles','on');
    h=findobj('Type','figure','Tag','D2-plot');
    set(0,'ShowHiddenHandles','off');  
    nr = size(h,1);
    set(fig,'KeyPressFcn','stop(''teststop'',gcbo)','Units','normalized');
    dat.nr=nr;dat.p2=[];dat.label=[];dat.traj=[];dat.labels=[];dat.trajs=[];
    set(fig,'color',[1 1 1],'UserData',dat);
    dath=get(h,'UserData');
    if nr~=1
        da=vertcat(dath{:,1}); datnr=sort(vertcat(da.nr));
        nr=find((datnr==(1:nr)')==0);
        if isempty(nr)
            nr=dat.nr;
        else
            dat.nr=nr(1); nr=nr(1);
            set(fig,'UserData',dat);
        end
    end
    MC.D2=[MC.D2 fig];
    nrfig=size(MC.D2,2);
    for siz=1:nrfig
        set(MC.D2(siz),'Name',sprintf('2Dplot:%d',siz));
    end
   	cla;
    %remember attributes?
    if (size(gds.plot2,1)<nr),drawing_attributes(gcf);end
    set(0,'showhiddenhandles','on')
    htmp = get(0, 'Children');
    set(0,'showhiddenhandles','off')
    % Generate a structure of handles to pass to callbacks, and store it.
    handles = guihandles(fig);
    guidata(fig, handles);
    if size(FigPos.D2,2)< nr, FigPos.D2{nr}=[0.5254 0.2826 0.4697 0.5677];end
    set(fig,'Position',FigPos.D2{nr});
    if (isfield(gds.plot2(nr,1),'lim')&& ~isempty(gds.plot2(nr,1).lim))
        if (gds.plot2(nr,1).lim(1,1))<(gds.plot2(nr,1).lim(1,2))
            set(gca,'XLIM',gds.plot2(nr,1).lim);
        else
            set(gca,'XLIM',[gds.plot2(nr,1).lim(1,2) gds.plot2(nr,1).lim(1,1)])
            set(gca,'XDir','reverse');
        end
    end
    if (isfield(gds.plot2(nr,2),'lim')&& ~isempty(gds.plot2(nr,2).lim))
        if (gds.plot2(nr,2).lim(1,1))<(gds.plot2(nr,2).lim(1,2))
            set(gca,'YLIM',gds.plot2(nr,2).lim);
        else
            set(gca,'YLIM',[gds.plot2(nr,2).lim(1,2) gds.plot2(nr,2).lim(1,1)])
            set(gca,'YDir','reverse');
        end
    end
    reset_callback(nrfig);
    if nargout > 0
		varargout{1} = fig;
    end
    figure(htmp(1));

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
		errordlg(lasterr,'error D2');
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
if (isfield(gds.plot2(nr,1),'lim')&& ~isempty(gds.plot2(nr,1).lim))
    if (gds.plot2(nr,1).lim(1,1))<(gds.plot2(nr,1).lim(1,2))
        set(gca,'XLIM',gds.plot2(nr,1).lim);
    else
        set(gca,'XLIM',[gds.plot2(nr,1).lim(1,2) gds.plot2(nr,1).lim(1,1)])
        set(gca,'XDir','reverse');
    end
end
if (isfield(gds.plot2(nr,2),'lim')&& ~isempty(gds.plot2(nr,2).lim))
    if (gds.plot2(nr,2).lim(1,1))<(gds.plot2(nr,2).lim(1,2))
        set(gca,'YLIM',gds.plot2(nr,2).lim);
    else
        set(gca,'YLIM',[gds.plot2(nr,2).lim(1,2) gds.plot2(nr,2).lim(1,1)])
        set(gca,'YDir','reverse');
    end
end
reset_callback(nr);

%----------------------------------------------------------------------
function redrawcurve(varargin)
global gds path_sys MC ; %lds;
if isnumeric(varargin{1}),nr=varargin{1};figure(MC.D2(nr));
else d=get(gcf,'UserData');nr=d.nr;figure(gcf);end
file=fullfile(path_sys,gds.system);save(file,'gds');
file=strcat(gds.curve.new,'.mat');
file=fullfile(path_sys,gds.system,gds.diagram,file);
if exist(file,'file')
    load(file);
    string=feval(gds.gui_load_draw);
    if isempty(strmatch(char(gds.plot2(nr,1).type),string,'exact'))|isempty(strmatch(char(gds.plot2(nr,2).type),string,'exact'))
        plot2=[inf inf];
    else
        plot2 = feval(gds.gui_make_ready,gds.plot2(nr,:),x);
    end
    feval(gds.gui_draw,2,plot2,file,'redraw');
else
    return;
end
%------------------------------------------------
function reset_callback(varargin)
global gds MC;
if (nargin==1)&&(isnumeric(varargin{1})),nr=varargin{1};h=get(MC.D2(nr),'CurrentAxes');
else d = get(gcf,'Name');nr=str2double(strrep(d,'2Dplot:',''));h=get(gcf,'CurrentAxes');end
axes(h);d=axis;
p=get(h,'View');
xmode = get(h,'Xdir');
ymode = get(h,'YDir');
cla reset;
user=get(MC.D2(nr),'UserData');
set(h,'View',p,'XDir',xmode,'YDir',ymode);axis(d);
if ~isfield(gds.plot2(user.nr,1),'label')
    gds.plot2(user.nr,1).label = [];
end
if ~isfield(gds.plot2(user.nr,2),'label')
    gds.plot2(user.nr,2).label = [];
end
xlabel(gds.plot2(user.nr,1).label);ylabel(gds.plot2(user.nr,2).label);

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
        plot2=feval(gds.gui_make_ready,gds.plot2(nr,:),x);
        feval(gds.gui_draw,2,plot2,filen,'redraw');
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
gds.type = ctype; gds.point = point;feval(strcat('gui_',deblank(ctype)));
gds.curve.new = start_filen;gds.curve.old = old_filen; 
file = fullfile(path_sys,gds.system); save(file,'gds');
feval(gds.gui_load_point,index,x,'not_save',start_file,1);


%----------------------------------------------------------------------
function close2D(varargin)
global gds path_sys MC calculation_progress FigPos nr
if calculation_progress~=0 
   return;
end;

d = find(MC.D2==gcf);
user = get(gcf,'UserData');
FigPos.D2{user.nr} = get(gcf,'Position');
h = findobj('Type','figure','Name',sprintf('2Dplot:%d',d));
delete(h);
h = findobj('Type','figure','Name',sprintf('Plotting region:%d',d));
delete(h);
MC.D2(d)=[];
if isempty(MC.D2)
%     gds.plot2 = [];
else
    gds.plot2(d,:) = [];
end
nrfig=size(MC.D2,2);
for siz=1:nrfig
    figure(MC.D2(siz));
    handles = guihandles(gcf);
    user = get(gcf,'UserData');
    nrf=get(MC.D2(siz),'Name');
    nr=str2num(strrep(nrf,'2Dplot:',''));
    c = findobj('Type','figure','Name',sprintf('Plotting region:%d',nr));
    user.nr = siz;
    set(MC.D2(siz),'Name',sprintf('2Dplot:%d',siz));
    set(c,'Name',sprintf('Plotting region:%d',siz));
    set(gcf,'UserData',user);
    guidata(gcf, handles);
end
