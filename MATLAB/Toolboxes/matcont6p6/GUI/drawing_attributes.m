function varargout = drawing_attributes(varargin)
% DRAWING_ATTRIBUTES Application M-file for drawing_attributes.fig
%    FIG = DRAWING_ATTRIBUTES launch drawing_attributes GUI.
%    DRAWING_ATTRIBUTES('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 08-Apr-2002 10:17:31
global gds oldplot nr driver_window;

if nargin == 1||nargin==0  % LAUNCH GUI
    if nargin==1,d=get(varargin{1},'UserData');nr = d.nr;
    else d=get(gcf,'UserData');nr=d.nr;end  
    fig = openfig(mfilename,'reuse');
	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
    oldplot=gds.plot2;
    load_draw(handles);  
	% Wait for callbacks to run and window to be dismissed:
% 	uiwait(fig);

	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK
	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
		errordlg(lasterr,mfilename);
        delete(driver_window);
	end

end


%--------------------------------------------------------------------
function list_Callback(h,eventdata,handles,varargin)
global gds nr;
switch eventdata
case 1
    h  = handles.abscissa;   h1 = handles.list_abscissa; h2 = handles.opt_abscissa; h3 = handles.eigs_abscissa;
case 2
    h  = handles.ordinate;    h1 = handles.list_ordinate; h2 = handles.opt_ordinate;  h3 = handles.eigs_ordinate;
end
index_selected = get(h1,'Value');
list_names = get(h1,'String');
name = list_names{index_selected};
index2 = get(h2,'Value');
index3 = get(h3,'Value');
list2 = get(h2,'String');
list3 = get(h3,'String');
switch eventdata
case 1
    gds.plot2(nr,1).type = name;
case 2
    gds.plot2(nr,2).type = name;
end
switch name
case 'Coordinates'
%     string = cellstr(char(gds.coordinates{:,1}));
    switch eventdata
        case 1
            set(handles.eigs_abscissa,'Visible','off');
        case 2
            set(handles.eigs_ordinate,'Visible','off');
    end
    [xxxx,plotopts] = feval(gds.gui_load_draw);
    if ~isempty(plotopts)
        switch eventdata
            case 1
                set(handles.opt_abscissa,'Visible','on');
            case 2
                set(handles.opt_ordinate,'Visible','on');
        end
            name2 = list2{index2};
            switch name2
                case 'Standard'
                    string=cellstr(char(gds.coordinates{:,1}));
                case 'Maximum'
                  tmp2 = char(gds.coordinates{:,1});
                  tmp1 = [];
                  tmp3 = [];
                  for i=1:size(gds.coordinates,1)
                      tmp1 = [tmp1; 'Max('];
                      tmp3 = [tmp3; ')'];
                  end
                  string = cellstr([tmp1 tmp2 tmp3]);
                case 'Minimum'
                  tmp2 = char(gds.coordinates{:,1});
                  tmp1 = [];
                  tmp3 = [];
                  for i=1:size(gds.coordinates,1)
                      tmp1 = [tmp1; 'Min('];
                      tmp3 = [tmp3; ')'];
                  end
                  string = cellstr([tmp1 tmp2 tmp3]);
                case 'Norm'
                  tmp2 = char(gds.coordinates{:,1});
                  tmp1 = [];
                  tmp3 = [];
                  for i=1:size(gds.coordinates,1)
                      tmp1 = [tmp1; 'Norm('];
                      tmp3 = [tmp3; ')'];
                  end
                  string = cellstr([tmp1 tmp2 tmp3]);
            end
    else
                    string=cellstr(char(gds.coordinates{:,1}));
    end
    set(h,'String',string,'Value',1);
case 'Parameters'
    string = cellstr(char(gds.parameters{:,1}));
    switch eventdata
        case 1
            set(handles.opt_abscissa,'Visible','off');
            set(handles.eigs_abscissa,'Visible','off');
        case 2
            set(handles.opt_ordinate,'Visible','off');
            set(handles.eigs_ordinate,'Visible','off');
    end
    set(h,'String',string,'Value',1);
case 'Stable Parameters'
    string = cellstr(char(gds.SParams{:,1}));
    switch eventdata
        case 1
            set(handles.opt_abscissa,'Visible','off');
            set(handles.eigs_abscissa,'Visible','off');
        case 2
            set(handles.opt_ordinate,'Visible','off');
            set(handles.eigs_ordinate,'Visible','off');
    end
    set(h,'String',string,'Value',1);
case 'Unstable Parameters'
    string = cellstr(char(gds.UParams{:,1}));
    switch eventdata
        case 1
            set(handles.opt_abscissa,'Visible','off');
            set(handles.eigs_abscissa,'Visible','off');
        case 2
            set(handles.opt_ordinate,'Visible','off');
            set(handles.eigs_ordinate,'Visible','off');
    end
    set(h,'String',string,'Value',1);
case 'Stable Parameter'
    str = gds.SParam{1,1};
    switch eventdata
        case 1
            set(handles.opt_abscissa,'Visible','off');
            set(handles.eigs_abscissa,'Visible','off');
        case 2
            set(handles.opt_ordinate,'Visible','off');
            set(handles.eigs_ordinate,'Visible','off');
    end
    set(h,'String',str,'Value',1);          
case 'Time'
      str = gds.time{1,1};
    switch eventdata
        case 1
            set(handles.opt_abscissa,'Visible','off');
            set(handles.eigs_abscissa,'Visible','off');
        case 2
            set(handles.opt_ordinate,'Visible','off');
            set(handles.eigs_ordinate,'Visible','off');
    end
      set(h,'String',str,'Value',1);
case 'Period'
    switch eventdata
        case 1
            set(handles.opt_abscissa,'Visible','off');
            set(handles.eigs_abscissa,'Visible','off');
        case 2
            set(handles.opt_ordinate,'Visible','off');
            set(handles.eigs_ordinate,'Visible','off');
    end
    set(h,'String',{'Period'},'Value',1);
case 'Multiplier'
    switch eventdata
        case 1
            set(handles.opt_abscissa,'Visible','off');
            set(handles.eigs_abscissa,'Visible','off');
        case 2
            set(handles.opt_ordinate,'Visible','off');
            set(handles.eigs_ordinate,'Visible','off');
    end
     str = feval(gds.gui_load_attributes);
     set(h,'String',str,'Value',1);
case 'Eigenvalues'
    switch eventdata
        case 1
            set(handles.opt_abscissa,'Visible','off');
            set(handles.eigs_abscissa,'Visible','on');
        case 2
            set(handles.opt_ordinate,'Visible','off');
            set(handles.eigs_ordinate,'Visible','on');
    end
    try
        name3 = list3{index3};
    catch
        name3 = 'Real part';
    end
    str = feval(gds.gui_load_attributes);
    switch name3
        case 'Real part'
            tmp1 = [];
            tmp3 = [];
            for i=1:gds.dim+1
                tmp1 = [tmp1; 'Re['];
                tmp3 = [tmp3; ']'];
            end
            str = cellstr([tmp1 char(str) tmp3]);
        case 'Imag. part'
            tmp1 = [];
            tmp3 = [];
            for i=1:gds.dim+1
                tmp1 = [tmp1; 'Im['];
                tmp3 = [tmp3; ']'];
            end
            str = cellstr([tmp1 char(str) tmp3]);
        case 'Norm'
            tmp1 = [];
            tmp3 = [];
            for i=1:gds.dim+1
                tmp1 = [tmp1; 'Norm['];
                tmp3 = [tmp3; ']'];
            end
            str = cellstr([tmp1 char(str) tmp3]);
     end
     set(h,'String',str,'Value',1);
case 'Userfunction'
    if gds.options.Userfunctions == 0
        error('No userfunctions are defined. Please choose another plotting option.');
    end
    switch eventdata
        case 1
            set(handles.opt_abscissa,'Visible','off');
            set(handles.eigs_abscissa,'Visible','off');
        case 2
            set(handles.opt_ordinate,'Visible','off');
            set(handles.eigs_ordinate,'Visible','off');
    end
    str = [];
    for i=1:length(gds.options.UserfunctionsInfo)
        str = [str; gds.options.UserfunctionsInfo(i).label];
    end
    str = cellstr(str);
     set(h,'String',str,'Value',1);
otherwise
     errordlg('Something went wrong');
end

%--------------------------------------------------------------------
function plotopt_Callback(h,eventdata,handles,varargin)
global gds nr;
switch eventdata
case 1
    h  = handles.abscissa;   h1 = handles.list_abscissa;    h2 = handles.opt_abscissa;
case 2
    h  = handles.ordinate;   h1 = handles.list_ordinate;    h2 = handles.opt_ordinate;
end
index_selected = get(h2,'Value');
list_names = get(h2,'String');
name = list_names{index_selected};
switch eventdata
case 1
    gds.plot2(nr,1).plotopt = name;
case 2
    gds.plot2(nr,2).plotopt = name;
end
switch name
case 'Standard'
    string = cellstr(char(gds.coordinates{:,1}));
    set(h,'String',string,'Value',1);
case 'Maximum'
      tmp2 = char(gds.coordinates{:,1});
      tmp1 = [];
      tmp3 = [];
      for i=1:size(gds.coordinates,1)
          tmp1 = [tmp1; 'Max('];
          tmp3 = [tmp3; ')'];
      end
      str = cellstr([tmp1 tmp2 tmp3]);
      set(h,'String',str,'Value',1);
case 'Minimum'
      tmp2 = char(gds.coordinates{:,1});
      tmp1 = [];
      tmp3 = [];
      for i=1:size(gds.coordinates,1)
          tmp1 = [tmp1; 'Min('];
          tmp3 = [tmp3; ')'];
      end
      str = cellstr([tmp1 tmp2 tmp3]);
      set(h,'String',str,'Value',1);
      set(h,'String',str,'Value',1);
case 'Norm'
      tmp2 = char(gds.coordinates{:,1});
      tmp1 = [];
      tmp3 = [];
      for i=1:size(gds.coordinates,1)
          tmp1 = [tmp1; 'Norm('];
          tmp3 = [tmp3; ')'];
      end
      str = cellstr([tmp1 tmp2 tmp3]);
      set(h,'String',str,'Value',1);
      set(h,'String',str,'Value',1);
otherwise
     errordlg('Something went wrong');
end

%--------------------------------------------------------------------
function eigsopt_Callback(h,eventdata,handles,varargin)
global gds nr;
switch eventdata
case 1
    h  = handles.abscissa;   h1 = handles.list_abscissa; h2 = handles.opt_abscissa; h3 = handles.eigs_abscissa;
case 2
    h  = handles.ordinate;    h1 = handles.list_ordinate; h2 = handles.opt_ordinate;  h3 = handles.eigs_ordinate;
end
index_selected = get(h3,'Value');
list_names = get(h3,'String');
name = list_names{index_selected};
switch eventdata
case 1
    gds.plot2(nr,1).eigsopt = name;
case 2
    gds.plot2(nr,2).eigsopt = name;
end
str = feval(gds.gui_load_attributes);
switch name
case 'Real part'
            tmp1 = [];
            tmp3 = [];
            for i=1:gds.dim+1
                tmp1 = [tmp1; 'Re['];
                tmp3 = [tmp3; ']'];
            end
            str = cellstr([tmp1 char(str) tmp3]);
      set(h,'String',str,'Value',1);
case 'Imag. part'
            tmp1 = [];
            tmp3 = [];
            for i=1:gds.dim+1
                tmp1 = [tmp1; 'Im['];
                tmp3 = [tmp3; ']'];
            end
            str = cellstr([tmp1 char(str) tmp3]);
      set(h,'String',str,'Value',1);
case 'Norm'
            tmp1 = [];
            tmp3 = [];
            for i=1:gds.dim+1
                tmp1 = [tmp1; 'Norm['];
                tmp3 = [tmp3; ']'];
            end
            str = cellstr([tmp1 char(str) tmp3]);
      set(h,'String',str,'Value',1);
otherwise
     errordlg('Something went wrong');
end
% --------------------------------------------------------------------
function ok_draw_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.ok_draw.
global gds nr path_sys;
gds.plot2(nr,1).val = get(handles.abscissa,'Value');
gds.plot2(nr,2).val = get(handles.ordinate,'Value');
xlabels = cellstr(get(handles.abscissa,'String'));
gds.plot2(nr,1).label = xlabels{gds.plot2(nr,1).val};
ylabels = cellstr(get(handles.ordinate,'String'));
gds.plot2(nr,2).label = ylabels{gds.plot2(nr,2).val};
file = fullfile(path_sys,gds.system);
save(file,'gds');
D2('OK',nr);
delete(handles.drawing_attributes);



%---------------------------------------------------------------------
function load_draw(handles)
global gds nr;
[string,plotopts] = feval(gds.gui_load_draw);
eigsopts = {'Real part'; 'Imag. part'; 'Norm'};
if (size(gds.plot2,1) >= nr)&&(size(gds.plot2,2)==2)&&~isempty(strmatch(char(gds.plot2(nr,1).type),string,'exact'))...
        &&~isempty(strmatch(char(gds.plot2(nr,2).type),string,'exact')) && isfield(gds.plot2(nr,1),'val')
    for d = 1:2
        str=[];
        switch d
        case 1
            h  = handles.abscissa;
            h1 = handles.list_abscissa;
            h2 = handles.opt_abscissa;
            h3 = handles.eigs_abscissa;
            set(handles.opt_abscissa,'Visible','off');
            set(handles.eigs_abscissa,'Visible','off');
        case 2
            h  = handles.ordinate;
            h1 = handles.list_ordinate;
            h2 = handles.opt_ordinate;
            h3 = handles.eigs_ordinate;
            set(handles.opt_ordinate,'Visible','off');
            set(handles.eigs_ordinate,'Visible','off');
        end
        name = gds.plot2(nr,d).type;
        if ~isfield(gds.plot2(nr,d),'plotopt') || isempty(gds.plot2(nr,d).plotopt)
            gds.plot2(nr,d).plotopt = 'Standard';
        end
        if ~isfield(gds.plot2(nr,d),'eigsopt') || isempty(gds.plot2(nr,d).eigsopt)
            gds.plot2(nr,d).eigsopt = 'Real part';
        end
        name2 = gds.plot2(nr,d).plotopt;
        name3 = gds.plot2(nr,d).eigsopt;
        if isempty(name3)
            name3 = 'Real part';
            gds.plot2(nr,d).eigsopt = 'Real part';
        end
        for j=1:size(string,1)
            if strcmp(name,string{j,1})
                break;
            end
        end
        for l=1:size(plotopts,1)
            if strcmp(name2,plotopts{l,1})
                break;
            end
        end
        for m=1:size(eigsopts,1)
            if strcmp(name3,eigsopts{m,1})
                break;
            end
        end
        switch name
        case 'Coordinates'
            switch d
                case 1
                    set(handles.opt_abscissa,'Visible','on');
                case 2
                    set(handles.opt_ordinate,'Visible','on');
            end
            switch name2
                case 'Standard'
                    str=cellstr(char(gds.coordinates{:,1}));
                case 'Maximum'
                  tmp2 = char(gds.coordinates{:,1});
                  tmp1 = [];
                  tmp3 = [];
                  for i=1:size(gds.coordinates,1)
                      tmp1 = [tmp1; 'Max('];
                      tmp3 = [tmp3; ')'];
                  end
                  str = cellstr([tmp1 tmp2 tmp3]);
                case 'Minimum'
                  tmp2 = char(gds.coordinates{:,1});
                  tmp1 = [];
                  tmp3 = [];
                  for i=1:size(gds.coordinates,1)
                      tmp1 = [tmp1; 'Min('];
                      tmp3 = [tmp3; ')'];
                  end
                  str = cellstr([tmp1 tmp2 tmp3]);
                case 'Norm'
                  tmp2 = char(gds.coordinates{:,1});
                  tmp1 = [];
                  tmp3 = [];
                  for i=1:size(gds.coordinates,1)
                      tmp1 = [tmp1; 'Norm('];
                      tmp3 = [tmp3; ')'];
                  end
                  str = cellstr([tmp1 tmp2 tmp3]);
            end
            if isempty(plotopts)
                str=cellstr(char(gds.coordinates{:,1}));
            end
            set(h2,'String',plotopts,'Value',l);
            set(h1,'String',string,'Value',j);
            set(h,'String',str,'Value',gds.plot2(nr,d).val);
        case 'Parameters'
            str = cellstr(char(gds.parameters{:,1}));
            set(h2,'String',plotopts,'Value',l);
            set(h1,'String',string,'Value',j);
            set(h,'String',str,'Value',gds.plot2(nr,d).val);
        case 'Stable Parameters'
            str = cellstr(char(gds.SParams{:,1}));
            set(h2,'String',plotopts,'Value',l);
            set(h1,'String',string,'Value',j);
            set(h,'String',str,'Value',gds.plot2(nr,d).val);
        case 'Unstable Parameters'
            str = cellstr(char(gds.UParams{:,1}));
            set(h2,'String',plotopts,'Value',l);
            set(h1,'String',string,'Value',j);
            set(h,'String',str,'Value',gds.plot2(nr,d).val);
        case 'Stable Parameter'
            str = gds.SParam{1,1};
            set(h2,'String',plotopts,'Value',l);
            set(h1,'String',string,'Value',j);
            set(h,'String',str,'Value',1);            
        case 'Time'
            str = gds.time{1,1};
            set(h2,'String',plotopts,'Value',l);
            set(h1,'String',string,'Value',j);
            set(h,'String',str,'Value',1);
        case 'Period'
            set(h2,'String',plotopts,'Value',l);
            set(h1,'String',string,'Value',j);
            set(h,'String',{'Period'},'Value',1);
        case 'eps1'
            set(h2,'String',plotopts,'Value',l);
            set(h1,'String',string,'Value',j);
            set(h,'String',{'eps1'},'Value',1);
            
        case 'Multiplier'
            str = feval(gds.gui_load_attributes);
            set(h2,'String',plotopts,'Value',l);
            set(h1,'String',string,'Value',j);
            set(h,'String',str,'Value',gds.plot2(nr,d).val);
        case 'Eigenvalues'
            switch d
                case 1
                    set(handles.eigs_abscissa,'Visible','on');
                case 2
                    set(handles.eigs_ordinate,'Visible','on');
            end
            str = feval(gds.gui_load_attributes);
            switch name3
            case 'Real part'
                        tmp1 = [];
                        tmp3 = [];
                        for i=1:gds.dim+1
                            tmp1 = [tmp1; 'Re['];
                            tmp3 = [tmp3; ']'];
                        end
                        str = cellstr([tmp1 char(str) tmp3]);
            case 'Imag. part'
                        tmp1 = [];
                        tmp3 = [];
                        for i=1:gds.dim+1
                            tmp1 = [tmp1; 'Im['];
                            tmp3 = [tmp3; ']'];
                        end
                        str = cellstr([tmp1 char(str) tmp3]);
            case 'Norm'
                        tmp1 = [];
                        tmp3 = [];
                        for i=1:gds.dim+1
                            tmp1 = [tmp1; 'Norm['];
                            tmp3 = [tmp3; ']'];
                        end
                        str = cellstr([tmp1 char(str) tmp3]);
            end
            set(h3,'String',eigsopts,'Value',m);
            set(h2,'String',plotopts,'Value',l);
            set(h1,'String',string,'Value',j);
            set(h,'String',str,'Value',gds.plot2(nr,d).val);  
case 'Userfunction'
    if gds.options.Userfunctions == 0
        error('No userfunctions are defined. Please choose another plotting option.');
    end
    str = [];
    for i=1:length(gds.options.UserfunctionsInfo)
        str = [str; gds.options.UserfunctionsInfo(i).label];
    end
    str = cellstr(str);
            set(h2,'String',plotopts,'Value',l);
            set(h1,'String',string,'Value',j);
     set(h,'String',str,'Value',gds.plot2(nr,d).val);  
        otherwise
            errordlg('Something went wrong');
        end
    end
else
    set(handles.list_abscissa,'String',string,'Value',1);
    set(handles.list_ordinate,'String',string,'Value',1);    
    str=cellstr(char(gds.coordinates{:,1}));
    set(handles.abscissa,'String',str,'Value',1);
    set(handles.ordinate,'String',str,'Value',1);
    set(handles.opt_ordinate,'String',plotopts,'Value',1);
    set(handles.eigs_ordinate,'String',eigsopts,'Value',1);
    set(handles.opt_abscissa,'String',plotopts,'Value',1);
    set(handles.eigs_abscissa,'String',eigsopts,'Value',1);
    gds.plot2(nr,1).type = 'Coordinates';
    gds.plot2(nr,2).type = 'Coordinates';
    gds.plot2(nr,1).plotopt = 'Standard';
    gds.plot2(nr,2).plotopt = 'Standard';
    set(handles.opt_abscissa,'Visible','on');
    set(handles.opt_ordinate,'Visible','on');
    set(handles.eigs_abscissa,'Visible','off');
    set(handles.eigs_ordinate,'Visible','off');
end
% --------------------------------------------------------------------
function cancel_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.cancel.
global gds oldplot nr;
if size(oldplot,1) < nr
    str=sprintf('Cancel isn''t possible\n Select the appropriate values and press the OK button');
    warndlg(str,'Cancel isn''t possible');
    return;
end
gds.plot2 = oldplot;
delete(handles.drawing_attributes);

