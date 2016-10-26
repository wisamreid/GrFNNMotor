function varargout= draw3_attributes(varargin)
% DRAW3_ATTRIBUTES Application M-file for draw3_attributes.fig
%    FIG = DRAW3_ATTRIBUTES launch draw3_attributes GUI.
%    DRAW3_ATTRIBUTES('callback_name', ...) invoke the named callback.

% Last Modified by GUIDE v2.0 08-Apr-2002 10:25:23
global gds oldplot nr driver_window;

if nargin == 0||nargin == 1  % LAUNCH GUI
    if nargin==1,d=get(varargin{1},'UserData');nr=d.nr;
    else d=get(gcf,'UserData');nr=d.nr;end  
	fig = openfig(mfilename,'reuse');
    % Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
    oldplot=gds.plot3;
    load_draw(handles);
	% Wait for callbacks to run and window to be dismissed:
	uiwait(fig);

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
    h=handles.ox;h1=handles.list_ox; h2 = handles.opt_ox; h3 = handles.eigs_ox;
case 2
    h=handles.oy;h1=handles.list_oy; h2 = handles.opt_oy; h3 = handles.eigs_oy;
case 3
    h=handles.oz;h1=handles.list_oz; h2 = handles.opt_oz; h3 = handles.eigs_oz;
end
index_selected = get(h1,'Value');
list_names= get(h1,'String');
name = list_names{index_selected};
index2 = get(h2,'Value');
index3 = get(h3,'Value');
list2 = get(h2,'String');
list3 = get(h3,'String');
switch eventdata
case 1
    gds.plot3(nr,1).type=name;
case 2
    gds.plot3(nr,2).type=name;
case 3
    gds.plot3(nr,3).type=name;
end
switch name
case 'Coordinates'
    switch eventdata
        case 1
            set(handles.eigs_ox,'Visible','off');
        case 2
            set(handles.eigs_oy,'Visible','off');
        case 3
            set(handles.eigs_oz,'Visible','off');
    end
    [xxxx,plotopts] = feval(gds.gui_load_draw);
    if ~isempty(plotopts)
        switch eventdata
        case 1
            set(handles.opt_ox,'Visible','on');
        case 2
            set(handles.opt_oy,'Visible','on');
        case 3
            set(handles.opt_oz,'Visible','on');
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
            set(handles.opt_ox,'Visible','off');
            set(handles.eigs_ox,'Visible','off');
        case 2
            set(handles.opt_oy,'Visible','off');
            set(handles.eigs_oy,'Visible','off');
        case 3
            set(handles.opt_oz,'Visible','off');
            set(handles.eigs_oz,'Visible','off');
    end
     set(h,'String',string,'value',1);
case 'Time'
      str=gds.time{1,1};
    switch eventdata
        case 1
            set(handles.opt_ox,'Visible','off');
            set(handles.eigs_ox,'Visible','off');
        case 2
            set(handles.opt_oy,'Visible','off');
            set(handles.eigs_oy,'Visible','off');
        case 3
            set(handles.opt_oz,'Visible','off');
            set(handles.eigs_oz,'Visible','off');
    end
      set(h,'String',str,'Value',1);
case 'Period'
    switch eventdata
        case 1
            set(handles.opt_ox,'Visible','off');
            set(handles.eigs_ox,'Visible','off');
        case 2
            set(handles.opt_oy,'Visible','off');
            set(handles.eigs_oy,'Visible','off');
        case 3
            set(handles.opt_oz,'Visible','off');
            set(handles.eigs_oz,'Visible','off');
    end
     set(h,'String',{'Period'},'Value',1);
case 'Multiplier'
    switch eventdata
        case 1
            set(handles.opt_ox,'Visible','off');
            set(handles.eigs_ox,'Visible','off');
        case 2
            set(handles.opt_oy,'Visible','off');
            set(handles.eigs_oy,'Visible','off');
        case 3
            set(handles.opt_oz,'Visible','off');
            set(handles.eigs_oz,'Visible','off');
    end
     str=feval(gds.gui_load_attributes);
     set(h,'String',str,'Value',1);
case 'Eigenvalues'
    switch eventdata
        case 1
            set(handles.opt_ox,'Visible','off');
            set(handles.eigs_ox,'Visible','on');
        case 2
            set(handles.opt_oy,'Visible','off');
            set(handles.eigs_oy,'Visible','on');
        case 3
            set(handles.opt_oz,'Visible','off');
            set(handles.eigs_oz,'Visible','on');
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
            set(handles.opt_ox,'Visible','off');
            set(handles.eigs_ox,'Visible','off');
        case 2
            set(handles.opt_oy,'Visible','off');
            set(handles.eigs_oy,'Visible','off');
        case 3
            set(handles.opt_oz,'Visible','off');
            set(handles.eigs_oz,'Visible','off');
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
    h=handles.ox;h1=handles.list_ox; h2 = handles.opt_ox; h3 = handles.eigs_ox;
case 2
    h=handles.oy;h1=handles.list_oy; h2 = handles.opt_oy; h3 = handles.eigs_oy;
case 3
    h=handles.oz;h1=handles.list_oz; h2 = handles.opt_oz; h3 = handles.eigs_oz;
end
index_selected = get(h2,'Value');
list_names = get(h2,'String');
name = list_names{index_selected};
switch eventdata
case 1
    gds.plot3(nr,1).plotopt = name;
case 2
    gds.plot3(nr,2).plotopt = name;
case 3
    gds.plot3(nr,3).plotopt = name;
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
    h=handles.ox;h1=handles.list_ox; h2 = handles.opt_ox; h3 = handles.eigs_ox;
case 2
    h=handles.oy;h1=handles.list_oy; h2 = handles.opt_oy; h3 = handles.eigs_oy;
case 3
    h=handles.oz;h1=handles.list_oz; h2 = handles.opt_oz; h3 = handles.eigs_oz;
end
index_selected = get(h3,'Value');
list_names = get(h3,'String');
name = list_names{index_selected};
switch eventdata
case 1
    gds.plot3(nr,1).eigsopt = name;
case 2
    gds.plot3(nr,2).eigsopt = name;
case 3
    gds.plot3(nr,3).eigsopt = name;
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
function ok_draw3_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.ok_draw.
global gds nr path_sys;
gds.plot3(nr,1).val=get(handles.ox,'Value');
gds.plot3(nr,2).val=get(handles.oy,'Value');
gds.plot3(nr,3).val=get(handles.oz,'Value');
oxs=cellstr(get(handles.ox,'String'));
gds.plot3(nr,1).label=oxs{gds.plot3(nr,1).val};
oys=cellstr(get(handles.oy,'String'));
gds.plot3(nr,2).label=oys{gds.plot3(nr,2).val};
ozs=cellstr(get(handles.oz,'String'));
gds.plot3(nr,3).label=ozs{gds.plot3(nr,3).val};
file=fullfile(path_sys,gds.system);
save(file,'gds');
plotD3('OK',nr);
delete(handles.draw3_attributes);

%--------------------------------------------------------------------
function load_draw(handles)
global gds nr;
[string,plotopts] = feval(gds.gui_load_draw);
eigsopts = {'Real part'; 'Imag. part'; 'Norm'};
if (size(gds.plot3,1) >= nr)&&(size(gds.plot3,2)==3)&&~isempty(strmatch(char(gds.plot3(nr,1).type),string,'exact'))&&~isempty(strmatch(char(gds.plot3(nr,2).type),string,'exact'))&~isempty(strmatch(char(gds.plot3(nr,3).type),string,'exact'))
    for d = 1:3
        str=[];
        switch d
        case 1
            h=handles.ox;h1=handles.list_ox;h2 = handles.opt_ox; h3 = handles.eigs_ox;
        case 2
            h=handles.oy;h1=handles.list_oy;h2 = handles.opt_oy; h3 = handles.eigs_oy;
        case 3
            h=handles.oz;h1=handles.list_oz;h2 = handles.opt_oz; h3 = handles.eigs_oz;
        end
        name = gds.plot3(nr,d).type;
        if ~isfield(gds.plot3(nr,d),'plotopt') || isempty(gds.plot3(nr,d).plotopt)
            gds.plot3(nr,d).plotopt = 'Standard';
        end
        if ~isfield(gds.plot3(nr,d),'eigsopt') || isempty(gds.plot3(nr,d).eigsopt)
            gds.plot3(nr,d).eigsopt = 'Real part';
        end
        name2 = gds.plot3(nr,d).plotopt;
        name3 = gds.plot3(nr,d).eigsopt;
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
                    set(handles.opt_ox,'Visible','on');
                case 2
                    set(handles.opt_oy,'Visible','on');
                case 3
                    set(handles.opt_oz,'Visible','on');
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
            set(h,'String',str,'Value',gds.plot3(nr,d).val);
        case 'Parameters'
            k=size(gds.parameters);
            str=cellstr(char(gds.parameters{:,1}));
            set(h1,'String',string,'Value',j);
            set(h,'String',str,'Value',gds.plot3(nr,d).val);
        case 'Time'
            str=gds.time{1,1};
            set(h1,'String',string,'Value',j);
            set(h,'String',str,'Value',1);
        case 'Period'
            set(h1,'String',string,'Value',j);
            set(h,'String',{'Period'},'Value',1);
        case 'Multiplier'
            str=feval(gds.gui_load_attributes);
            set(h1,'String',string,'Value',j);
            set(h,'String',str,'Value',gds.plot3(nr,d).val);
        case 'Eigenvalues'
            switch d
                case 1
                    set(handles.eigs_ox,'Visible','on');
                case 2
                    set(handles.eigs_oy,'Visible','on');
                case 3
                    set(handles.eigs_oz,'Visible','on');
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
            set(h,'String',str,'Value',gds.plot3(nr,d).val);  
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
     set(h,'String',str,'Value',1);  
        otherwise
            errordlg('Something went wrong');
        end
    end
else
    set(handles.list_ox,'String',string,'Value',1);
    set(handles.list_oy,'String',string,'Value',1);
    set(handles.list_oz,'String',string,'Value',1);
    str=cellstr(char(gds.coordinates{:,1}));
    set(handles.ox,'String',str,'Value',1);
    set(handles.oy,'String',str,'Value',1);
    set(handles.oz,'String',str,'Value',1);
    gds.plot3(nr,1).type='Coordinates';
    gds.plot3(nr,2).type='Coordinates';
    gds.plot3(nr,3).type='Coordinates';
end




% --------------------------------------------------------------------
function cancel_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.cancel.
global gds oldplot nr;
if size(oldplot,1)<nr
    str=sprintf('Cancel isn''t possible\n Select the appropriate values and press the OK button');
    warndlg(str,'Cancel isn''t possible');
    return;
end
gds.plot3=oldplot;
delete(handles.draw3_attributes);
