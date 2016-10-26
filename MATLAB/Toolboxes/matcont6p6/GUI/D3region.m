function varargout = D3region(varargin)
% D3REGION Application M-file for D3region.fig
%    FIG = D3REGION launch D3region GUI.
%    D3REGION('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 02-Dec-2003 09:12:05
global gds path_sys MC calculation_progress FigPos 
if ~ischar(varargin{1})  % LAUNCH GUI
    d = get(varargin{1},'UserData');nr1=d.nr;
    d = get(varargin{1},'Name');
    nr=str2double(strrep(d,'3Dplot:',''));
    set(0,'ShowHiddenHandles','on');
    c = findobj('Type','figure','Name',sprintf('3D Plotting region:%d',nr));
    if ~isempty(c),figure(c);return;end
    set(0,'ShowHiddenHandles','off');
	fig = openfig(mfilename,'new');
	set(fig,'Color',[1 1 1]);
    set(fig,'Name',sprintf('3D Plotting region:%d',nr));

	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
    if (isfield(gds.plot3,'lim')&& ~isempty(gds.plot3(nr1,1).lim))
        set(handles.x1,'String',num2str(gds.plot3(nr1,1).lim(1,1)));
        set(handles.x2,'string',num2str(gds.plot3(nr1,1).lim(1,2)));
    end
    if (isfield(gds.plot3,'lim')&& ~isempty(gds.plot3(nr1,2).lim))
        set(handles.y1,'string',num2str(gds.plot3(nr1,2).lim(1,1)));
        set(handles.y2,'string',num2str(gds.plot3(nr1,2).lim(1,2)));
    end
    if (isfield(gds.plot3,'lim')&& ~isempty(gds.plot3(nr1,3).lim))
        set(handles.z1,'string',num2str(gds.plot3(nr1,3).lim(1,1)));
        set(handles.z2,'string',num2str(gds.plot3(nr1,3).lim(1,2)));
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
        figure(tmp);
	catch
		disp(lasterr);
	end

end



% --------------------------------------------------------------------
function x1_Callback(h, eventdata, handles, varargin)
global gds MC 
d = get(gcf,'Name');
nr=str2double(strrep(d,'3D Plotting region:',''));
h=get(MC.D3(nr),'CurrentAxes');
a=get(h,'XLIM');
mode = get(h,'XDir');
b = str2double(get(handles.x1,'String'));
if strcmp(mode,'normal')
    if (a(1,2)<b)
        set(h,'XDir','reverse');
        a(1,1)=b;
        set(h,'XLIM',a(end:-1:1));
        gds.plot3(nr,1).lim=a;
    elseif (a(1,2)==b)
        errordlg('Values must be increasing');
        set(handles.x1,'String',num2str(a(1,1)));
        return
    elseif (a(1,2)>b)
        a(1,1)=b;
        set(h,'XLIM',a);
        gds.plot3(nr,1).lim=a;        
    end
else
    if (a(1,1)> b)
        set(h,'XDir','normal');
        a(1,2) = b;
        set(h,'XLIM',a(end:-1:1));
        gds.plot3(nr,1).lim=a(end:-1:1);
    elseif (a(1,1)==b)
        errordlg('Values must be increasing');
        set(handles.x1,'String',num2str(a(1,1)));
        return
    elseif (a(1,1)<b)
        a(1,2)=b;
        set(h,'XLIM',a);
        gds.plot3(nr,1).lim=a(end:-1:1);        
    end
end 
plotD3('reset_callback',nr);
plotD3('redrawcurve',nr);




% --------------------------------------------------------------------
function y1_Callback(h, eventdata, handles, varargin)
global gds MC
d = get(gcf,'Name');
nr=str2double(strrep(d,'3D Plotting region:',''));
h=get(MC.D3(nr),'CurrentAxes');
a = get(h,'YLIM');
mode = get(h,'YDir');
b=str2double(get(handles.y1,'String'));
if strcmp(mode,'normal')
    if (a(1,2)<b)
        set(h,'YDir','reverse');
        a(1,1)=b;
        set(h,'YLIM',a(end:-1:1));
        gds.plot3(nr,2).lim=a;
    elseif (a(1,2)==b)
        errordlg('Values must be increasing');
        set(handles.y1,'String',num2str(a(1,1)));
        return
    elseif (a(1,2)>b)
        a(1,1)=b;
        set(h,'YLIM',a);
        gds.plot3(nr,2).lim=a;        
    end
else
    if (a(1,1)> b)
        set(h,'YDir','normal');
        a(1,2) = b;
        set(h,'YLIM',a(end:-1:1));
        gds.plot3(nr,2).lim=a(end:-1:1);
    elseif (a(1,1)==b)
        errordlg('Values must be increasing');
        set(handles.y1,'String',num2str(a(1,1)));
        return
    elseif (a(1,1)<b)
        a(1,2)=b;
        set(h,'YLIM',a);
        gds.plot3(nr,2).lim=a(end:-1:1);        
    end
end 
plotD3('reset_callback',nr);
plotD3('redrawcurve',nr);




% --------------------------------------------------------------------
function x2_Callback(h, eventdata, handles, varargin)
global gds MC
d = get(gcf,'Name');
nr=str2double(strrep(d,'3D Plotting region:',''));
h=get(MC.D3(nr),'CurrentAxes');
a=get(h,'XLIM');
mode = get(h,'XDir');
b=str2double(get(handles.x2,'String'));
if strcmp(mode,'normal')
    if (a(1,1)>b)
        set(h,'XDir','reverse');
        a(1,2)=b;
        set(h,'XLIM',a(end:-1:1));
        gds.plot3(nr,1).lim=a;
    elseif (a(1,2)==b)
        errordlg('Values must be increasing');
        set(handles.x2,'String',num2str(a(1,1)));
        return
    else
        a(1,2)=b;
        set(h,'XLIM',a);
        gds.plot3(nr,1).lim=a;        
    end
else
    if (a(1,2)< b)
        set(h,'XDir','normal');
        a(1,1) = b;
        set(h,'XLIM',a(end:-1:1));
        gds.plot3(nr,1).lim=a(end:-1:1);
    elseif (a(1,2)==b)
        errordlg('Values must be increasing');
        set(handles.x2,'String',num2str(a(1,1)));
        return
    else
        a(1,1)=b;
        set(h,'XLIM',a);
        gds.plot3(nr,1).lim=a(end:-1:1);        
    end
end 
plotD3('reset_callback',nr);
plotD3('redrawcurve',nr);




% --------------------------------------------------------------------
function y2_Callback(h, eventdata, handles, varargin)
global gds MC 
d = get(gcf,'Name');
nr=str2double(strrep(d,'3D Plotting region:',''));
h=get(MC.D3(nr),'CurrentAxes');
a = get(h,'YLIM');
mode = get(h,'YDir');
b=str2double(get(handles.y2,'String'));
if strcmp(mode,'normal')
    if (a(1,1)>b)
        set(h,'YDir','reverse');
        a(1,2)=b;
        set(h,'YLIM',a(end:-1:1));
        gds.plot3(nr,2).lim=a;
    elseif (a(1,2)==b)
        errordlg('Values must be increasing');
        set(handles.y2,'String',num2str(a(1,1)));
        return
    else
        a(1,2)=b;
        set(h,'YLIM',a);
        gds.plot3(nr,2).lim=a;        
    end
else
    if (a(1,2)< b)
        set(h,'YDir','normal');
        a(1,1) = b;
        set(h,'YLIM',a(end:-1:1));
        gds.plot3(nr,2).lim=a(end:-1:1);
    elseif (a(1,2)==b)
        errordlg('Values must be increasing');
        set(handles.y2,'String',num2str(a(1,1)));
        return
    else
        a(1,1)=b;
        set(h,'YLIM',a);
        gds.plot3(nr,2).lim=a(end:-1:1);        
    end
end 
plotD3('reset_callback',nr);
plotD3('redrawcurve',nr);
% --------------------------------------------------------------------
function z2_Callback(h, eventdata, handles, varargin)
global gds MC 
d = get(gcf,'Name');
nr=str2double(strrep(d,'3D Plotting region:',''));
h=get(MC.D3(nr),'CurrentAxes');
a=get(h,'ZLIM');
mode = get(h,'ZDir');
b=str2double(get(handles.z2,'String'));
if strcmp(mode,'normal')
    if (a(1,1)>b)
        set(h,'ZDir','reverse');
        a(1,2)=b;
        set(h,'ZLIM',a(end:-1:1));
        gds.plot3(nr,3).lim=a;
    elseif (a(1,2)==b)
        errordlg('Values must be increasing');
        set(handles.z2,'String',num2str(a(1,1)));
        return
    else
        a(1,2)=b;
        set(h,'ZLIM',a);
        gds.plot3(nr,3).lim=a;        
    end
else
    if (a(1,2)< b)
        set(h,'ZDir','normal');
        a(1,1) = b;
        set(h,'ZLIM',a(end:-1:1));
        gds.plot3(nr,3).lim=a(end:-1:1);
    elseif (a(1,2)==b)
        errordlg('Values must be increasing');
        set(handles.z2,'String',num2str(a(1,1)));
        return
    else
        a(1,1)=b;
        set(h,'ZLIM',a);
        gds.plot3(nr,3).lim=a(end:-1:1);        
    end
end 
plotD3('reset_callback',nr);
plotD3('redrawcurve',nr);

% --------------------------------------------------------------------
function z1_Callback(h, eventdata, handles, varargin)
global gds MC 
d = get(gcf,'Name');
nr=str2double(strrep(d,'3D Plotting region:',''));
h=get(MC.D3(nr),'CurrentAxes');
a = get(h,'ZLIM');
mode = get(h,'ZDir');
b=str2double(get(handles.z1,'String'));
if strcmp(mode,'normal')
    if (a(1,2)<b)
        set(h,'ZDir','reverse');
        a(1,1)=b;
        set(h,'ZLIM',a(end:-1:1));
        gds.plot3(nr,3).lim=a;
    elseif (a(1,2)==b)
        errordlg('Values must be increasing');
        set(handles.z1,'String',num2str(a(1,1)));
        return
    elseif (a(1,2)>b)
        a(1,1)=b;
        set(h,'ZLIM',a);
        gds.plot3(nr,3).lim=a;        
    end
else
    if (a(1,1)> b)
        set(h,'ZDir','normal');
        a(1,2) = b;
        set(h,'ZLIM',a(end:-1:1));
        gds.plot3(nr,3).lim=a(end:-1:1);
    elseif (a(1,1)==b)
        errordlg('Values must be increasing');
        set(handles.z1,'String',num2str(a(1,1)));
        return
    elseif (a(1,1)<b)
        a(1,2)=b;
        set(h,'ZLIM',a);
        gds.plot3(nr,3).lim=a(end:-1:1);        
    end
end 
plotD3('reset_callback',nr);
plotD3('redrawcurve',nr);

