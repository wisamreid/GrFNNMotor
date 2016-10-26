function varargout = D2region(varargin)
% D2REGION Application M-file for D2region.fig
%    FIG = D2REGION launch D2region GUI.
%    D2REGION('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 02-Dec-2003 08:31:25
global gds 
if ~ischar(varargin{1})  % LAUNCH GUI
    d = get(varargin{1},'UserData');    nr1 = d.nr;
    d = get(varargin{1},'Name');
    nr=str2double(strrep(d,'2Dplot:',''));
    set(0,'ShowHiddenHandles','on');
    c = findobj('Type','figure','Name',sprintf('2D Plotting region:%d',nr));
    if ~isempty(c),figure(c);return;end
    set(0,'ShowHiddenHandles','off');
	fig = openfig(mfilename,'new');
	% Use system color scheme for figure:
	set(fig,'Color',[1 1 1]);
    set(fig,'Name',sprintf('2D Plotting region:%d',nr));
    
	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
    if (isfield(gds.plot2(nr1,1),'lim')&& ~isempty(gds.plot2(nr1,1).lim))
        set(handles.abs1,'String',num2str(gds.plot2(nr1,1).lim(1,1)));
        set(handles.abs2,'string',num2str(gds.plot2(nr1,1).lim(1,2)));
    end
    if (isfield(gds.plot2(nr1,2),'lim')&& ~isempty(gds.plot2(nr1,2).lim))
        set(handles.ord1,'string',num2str(gds.plot2(nr1,2).lim(1,1)));
        set(handles.ord2,'string',num2str(gds.plot2(nr1,2).lim(1,2)));
    end

	guidata(fig, handles);

	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
        tmp = gcf;
        if (nargout)
			[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
		else
			feval(varargin{:}); % FEVAL switchyard
        end
        figure(tmp)
	catch
		disp(lasterr);
	end

end



% --------------------------------------------------------------------
function abs1_Callback(h, eventdata, handles, varargin)
global gds MC 
d = get(gcf,'Name');
nr=str2double(strrep(d,'2D Plotting region:',''));
h=get(MC.D2(nr),'CurrentAxes');
a=get(h,'XLIM');
mode=get(h,'Xdir');
b=str2double(get(handles.abs1,'String'));
if strcmp(mode,'normal')
    if (a(1,2)<b)
        set(h,'XDir','reverse');
        a(1,1)=b;
        set(h,'XLIM',a(end:-1:1));
        gds.plot2(nr,1).lim=a;
    elseif (a(1,2)==b)
        errordlg('Values must be increasing');
        set(handles.abs1,'String',num2str(a(1,1)));
        return
    elseif (a(1,2)>b)
        a(1,1)=b;
        set(h,'XLIM',a);
        gds.plot2(nr,1).lim=a;        
    end
else
    if (a(1,1)> b)
        set(h,'XDir','normal');
        a(1,2) = b;
        set(h,'XLIM',a(end:-1:1));
        gds.plot2(nr,1).lim=a(end:-1:1);
    elseif (a(1,1)==b)
        errordlg('Values must be increasing');
        set(handles.abs1,'String',num2str(a(1,1)));
        return
    elseif (a(1,1)<b)
        a(1,2)=b;
        set(h,'XLIM',a);
        gds.plot2(nr,1).lim=a(end:-1:1);        
    end
end 
D2('reset_callback',nr);
D2('redrawcurve',nr);

% --------------------------------------------------------------------
function ord1_Callback(h, eventdata, handles, varargin)

global gds MC
d = get(gcf,'Name');
nr=str2double(strrep(d,'2D Plotting region:',''));
h=get(MC.D2(nr),'CurrentAxes');
a=get(h,'YLIM');
b=str2double(get(handles.ord1,'String'));
mode=get(h,'YDir');
b=str2double(get(handles.ord1,'String'));
if strcmp(mode,'normal')
    if (a(1,2)<b)
        set(h,'YDir','reverse');
        a(1,1)=b;
        set(h,'YLIM',a(end:-1:1));
        gds.plot2(nr,2).lim=a;
    elseif (a(1,2)==b)
        errordlg('Values must be increasing');
        set(handles.ord1,'String',num2str(a(1,1)));
        return
    elseif (a(1,2)>b)
        a(1,1)=b;
        set(h,'YLIM',a);
        gds.plot2(nr,2).lim=a;        
    end
else
    if (a(1,1)> b)
        set(h,'YDir','normal');
        a(1,2) = b;
        set(h,'YLIM',a(end:-1:1));
        gds.plot2(nr,2).lim=a(end:-1:1);
    elseif (a(1,1)==b)
        errordlg('Values must be increasing');
        set(handles.ord1,'String',num2str(a(1,1)));
        return
    elseif (a(1,1)<b)
        a(1,2)=b;
        set(h,'YLIM',a);
        gds.plot2(nr,2).lim=a(end:-1:1);        
    end
end 
D2('reset_callback',nr);
D2('redrawcurve',nr);

% --------------------------------------------------------------------
function abs2_Callback(h, eventdata, handles, varargin)

global gds MC
d = get(gcf,'Name');
nr=str2double(strrep(d,'2D Plotting region:',''));
h=get(MC.D2(nr),'CurrentAxes');
a=get(h,'XLIM');
mode = get(h,'XDir');
b=str2double(get(handles.abs2,'String'));
if strcmp(mode,'normal')
    if (a(1,1)>b)
        set(h,'XDir','reverse');
        a(1,2)=b;
        set(h,'XLIM',a(end:-1:1));
        gds.plot2(nr,1).lim=a;
    elseif (a(1,2)==b)
        errordlg('Values must be increasing');
        set(handles.abs2,'String',num2str(a(1,1)));
        return
    else
        a(1,2)=b;
        set(h,'XLIM',a);
        gds.plot2(nr,1).lim=a;        
    end
else
    if (a(1,2)< b)
        set(h,'XDir','normal');
        a(1,1) = b;
        set(h,'XLIM',a(end:-1:1));
        gds.plot2(nr,1).lim=a(end:-1:1);
    elseif (a(1,2)==b)
        errordlg('Values must be increasing');
        set(handles.abs2,'String',num2str(a(1,1)));
        return
    else
        a(1,1)=b;
        set(h,'XLIM',a);
        gds.plot2(nr,1).lim=a(end:-1:1);        
    end
end 
D2('reset_callback',nr);
D2('redrawcurve',nr);

% --------------------------------------------------------------------
function ord2_Callback(h, eventdata, handles, varargin)

global gds MC 
d = get(gcf,'Name');
nr=str2double(strrep(d,'2D Plotting region:',''));
h=get(MC.D2(nr),'CurrentAxes');
a=get(h,'YLIM');
mode = get(h,'YDir');
b=str2double(get(handles.ord2,'String'));
if strcmp(mode,'normal')
    if (a(1,1)>b)
        set(h,'YDir','reverse');
        a(1,2)=b;
        set(h,'YLIM',a(end:-1:1));
        gds.plot2(nr,2).lim=a;
    elseif (a(1,2)==b)
        errordlg('Values must be increasing');
        set(handles.ord2,'String',num2str(a(1,1)));
        return
    else
        a(1,2)=b;
        set(h,'YLIM',a);
        gds.plot2(nr,2).lim=a;        
    end
else
    if (a(1,2)< b)
        set(h,'YDir','normal');
        a(1,1) = b;
        set(h,'YLIM',a(end:-1:1));
        gds.plot2(nr,2).lim=a(end:-1:1);
    elseif (a(1,2)==b)
        errordlg('Values must be increasing');
        set(handles.ord2,'String',num2str(a(1,1)));
        return
    else
        a(1,1)=b;
        set(h,'YLIM',a);
        gds.plot2(nr,2).lim=a(end:-1:1);        
    end
end 
D2('reset_callback',nr);
D2('redrawcurve',nr);