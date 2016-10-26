function varargout = plotD3(varargin)
% PLOTD3 Application M-file for plotD3.fig
%    FIG = PLOTD3 launch plotD3 GUI.
%    PLOTD3('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 13-Dec-2001 09:06:13
global gds MC driver_window FigPos;
if nargin == 0  % LAUNCH GUI
    fig = openfig(mfilename,'new');  
%     set(fig,'KeyPressFcn','stop(''teststop'',gcbo)','Units','normalized');
    set(0,'ShowHiddenHandles','on');
    h=findobj('Type','figure','Tag','D3-plot');
    set(0,'ShowHiddenHandles','off');    
    nr = size(h,1);    
    set(fig,'KeyPressFcn','stop(''teststop'',gcbo)','Units','normalized');
    dat.nr=nr;dat.p3=[];dat.label=[];dat.traj=[];dat.labels=[];dat.trajs=[];
    set(fig,'color',[1 1 1],'UserData',dat);
    dath=get(h,'UserData');
    if nr~=1
        da=vertcat(dath{:,1});datnr=sort(vertcat(da.nr));
        nr=find((datnr==(1:nr)')==0);
        if isempty(nr)
            nr=dat.nr;
        else
            dat.nr=nr(1);nr=nr(1);
            set(fig,'UserData',dat);
        end
    end
    MC.D3=[MC.D3 fig];
    nrfig=size(MC.D3,2);
    for siz=1:nrfig
        set(MC.D3(siz),'Name',sprintf('3Dplot:%d',siz));
    end
    cla;
    %remember attributes?
    if (size(gds.plot3,1)<nr),draw3_attributes(gcf);end
	% Generate a structure of handles to pass to callbacks, and store it.   
    handles = guihandles(fig);
    guidata(fig, handles);
    if size(FigPos.D3,2)< nr, FigPos.D3{nr}=[0.004 0.335 0.459 0.439];end
    set(fig,'Position',FigPos.D3{nr});
    if (isfield(gds.plot3,'lim')&& ~isempty(gds.plot3(nr,1).lim))
        if gds.plot3(nr,1).lim(1,1)>gds.plot3(nr,1).lim(1,2)
            set(gca,'XDir','reverse');
            set(gca,'XLIM',[gds.plot3(nr,1).lim(1,2) gds.plot3(nr,1).lim(1,1)]);
        else
            set(gca,'XLIM',gds.plot3(nr,1).lim);
        end
    end
    if (isfield(gds.plot3,'lim')&& ~isempty(gds.plot3(nr,2).lim))
        if gds.plot3(nr,2).lim(1,1)>gds.plot3(nr,2).lim(1,2)
            set(gca,'YDir','reverse');
            set(gca,'YLIM',[gds.plot3(nr,2).lim(1,2) gds.plot3(nr,2).lim(1,1)]);
        else
            set(gca,'YLIM',gds.plot3(nr,2).lim);
        end
    end
    if (isfield(gds.plot3,'lim')&& ~isempty(gds.plot3(nr,3).lim))
        if gds.plot3(nr,3).lim(1,1)>gds.plot3(nr,3).lim(1,2)
            set(gca,'ZDir','reverse');
            set(gca,'ZLIM',[gds.plot3(nr,3).lim(1,2) gds.plot3(nr,3).lim(1,1)]);
        else
            set(gca,'ZLIM',gds.plot3(nr,3).lim);
        end
    end
    reset_callback(nrfig);
	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
		errordlg(lasterr,'error D3');
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
if (isfield(gds.plot3,'lim')&& ~isempty(gds.plot3(nr,1).lim))
    if gds.plot3(nr,1).lim(1,1)>gds.plot3(nr,1).lim(1,2)
        set(gca,'XDir','reverse');
        set(gca,'XLIM',[gds.plot3(nr,1).lim(1,2) gds.plot3(nr,1).lim(1,1)]);
    else
        set(gca,'XLIM',gds.plot3(nr,1).lim);
    end
end
if (isfield(gds.plot3,'lim')&& ~isempty(gds.plot3(nr,2).lim))
    if gds.plot3(nr,2).lim(1,1)>gds.plot3(nr,2).lim(1,2)
        set(gca,'YDir','reverse');
        set(gca,'YLIM',[gds.plot3(nr,2).lim(1,2) gds.plot3(nr,2).lim(1,1)]);
    else
        set(gca,'YLIM',gds.plot3(nr,2).lim);
    end
end
if (isfield(gds.plot3,'lim')&& ~isempty(gds.plot3(nr,3).lim))
    if gds.plot3(nr,3).lim(1,1)>gds.plot3(nr,3).lim(1,2)
        set(gca,'ZDir','reverse');
        set(gca,'ZLIM',[gds.plot3(nr,3).lim(1,2) gds.plot3(nr,3).lim(1,1)]);
    else
        set(gca,'ZLIM',gds.plot3(nr,3).lim);
    end
end
reset_callback(nr); 

%-----------------------------------------------------------------
function redrawcurve(varargin)
global gds path_sys lds MC;
if isnumeric(varargin{1}),nr=varargin{1};figure(MC.D3(nr));
else d=get(gcf,'UserData');nr=d.nr;figure(gcf);end
file=fullfile(path_sys,gds.system);save(file,'gds');
file=strcat(gds.curve.new,'.mat');
file=fullfile(path_sys,gds.system,gds.diagram,file);
if exist(file,'file')
    load(file);
    feval(strcat('gui_',deblank(gds.type)));
    string=feval(gds.gui_load_draw);
    if isempty(strmatch(char(gds.plot3(nr,1).type),string,'exact'))|isempty(strmatch(char(gds.plot3(nr,2).type),string,'exact'))|isempty(strmatch(char(gds.plot3(nr,3).type),string,'exact'))        
        plo = [inf inf inf] ;   
    else
        plo = feval(gds.gui_make_ready,gds.plot3(nr,:),x);
    end
    feval(gds.gui_draw,3,plo,file,'redraw');
else
    return;
end


%------------------------------------------------
function reset_callback(varargin)
global gds MC;
if (nargin==1)&&(isnumeric(varargin{1})),nr=varargin{1};h=get(MC.D3(nr),'CurrentAxes');
else d=get(gcf,'Name');nr=str2double(strrep(d,'3Dplot:',''));h=get(gcf,'CurrentAxes');end
axes(h);d=axis;
p=get(h,'View');
xmode = get(h,'Xdir');
ymode = get(h,'YDir');
zmode = get(h,'ZDir');
cla reset;user=get(MC.D3(nr),'UserData');
set(h,'View',p,'XDir',xmode,'YDir',ymode,'ZDir',zmode);axis(d);
xlabel(gds.plot3(user.nr,1).label);ylabel(gds.plot3(user.nr,2).label);zlabel(gds.plot3(user.nr,3).label);
%----------------------------------------------------------------------
function redrawdiagram(varargin)
global gds path_sys %lds;
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
        gds.type = ctype; gds.point = point;feval(strcat('gui_',deblank(ctype)));
        gds.curve.new = strrep(w.mat{i},'.mat','');gds.curve.old = strrep(w.mat{i},'.mat',''); 
        file = fullfile(path_sys,gds.system);
        save(file,'gds');
        feval(gds.gui_load_point,index,x,'not_save',filen,1);      
        plot3=feval(gds.gui_make_ready,gds.plot3(nr,:),x);
        feval(gds.gui_draw,3,plot3,filen,'redraw');
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
function close3D(varargin)
global gds path_sys MC calculation_progress  FigPos;%lds;
if calculation_progress~=0 
   return;
end;
d = find(MC.D3==gcf);
user = get(MC.D3(d),'UserData');
FigPos.D3{user.nr} = get(MC.D3(d),'Position');
h = findobj('Type','figure','Name',sprintf('3D Plotting region:%d',user.nr));
delete(h);
h = findobj('Type','figure','Name',sprintf('3Dplot:%d',d));
delete(h);
MC.D3(d) = [];
if isempty(MC.D3)
%     gds.plot2 = [];
else
    gds.plot3(d,:) = [];
end
nrfig=size(MC.D3,2);
for siz=1:nrfig
    figure(MC.D3(siz));
    handles = guihandles(gcf);
    user = get(gcf,'UserData');
    nrf=get(MC.D3(siz),'Name');
    nr=str2num(strrep(nrf,'3Dplot:',''));
    c = findobj('Type','figure','Name',sprintf('3D Plotting region:%d',nr));
    user.nr = siz;
    set(MC.D3(siz),'Name',sprintf('3Dplot:%d',siz));
    set(c,'Name',sprintf('3D Plotting region:%d',siz));
    set(gcf,'UserData',user);
    guidata(gcf, handles);
end