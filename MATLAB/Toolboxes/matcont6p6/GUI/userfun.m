function varargout = userfun(varargin)
% USERFUN Application M-file for userfun.fig
%    FIG = USERFUN launch userfun GUI.
%    USERFUN('callback_name', ...) invoke the named callback.


% Last Modified by GUIDE v2.0 04-Sep-2002 11:20:31
global gds oldgds path_sys driver_window MC;
if nargin == 0  % LAUNCH GUI

	fig = openfig(mfilename,'reuse');
    oldgds = gds;
	% Generate a structure of handles to pass to callbacks, and store it. 
	handles = guihandles(fig);
	guidata(fig, handles);
    load_listbox1(handles);
    %starter-window open?
    if ~isempty(MC.starter),gds.open.figuur=1;else gds.open.figuur=0;end
    delete(MC.starter);MC.starter=[];
    % numeric window open?   
    if ~isempty(MC.numeric_fig), gds.open.numeric_fig=1;else gds.open.numeric_fig=0;end
    close(MC.numeric_fig);MC.numeric_fig=[];
	% Wait for callbacks to run and window to be dismissed:
	uiwait(fig);

	if nargout > 0
		varargout{1} = fig;
	end

elseif ischar(varargin{1}) % INVOKE NAMED SUBFUNCTION OR CALLBACK

	try
		[varargout{1:nargout}] = feval(varargin{:}); % FEVAL switchyard
	catch
		disp(lasterr);      
	end
end

% --------------------------------------------------------------------
function adbutton_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.adbutton.

global gds;
label = get(handles.label,'String');
name  = get(handles.name,'String');
if length(label)>2
    warndlg('A label exist of 2 characters')
    set(handles.label,'String',label(1:2));
elseif length(label)<2
    while length(label)<2
        label = sprintf('%s%c',label,char(32));
    end
    set(handles.label,'String',label(1:2));
end
if length(name)<1
    warndlg('You have to enter the name of the function');
    return;
end
if isfield(gds,'userfunction')%userfunctions already exists
    dim = size(gds.userfunction,2);
    namestr  = cellstr(char(gds.options.UserfunctionsInfo.name));
    labelstr = cellstr(char(gds.options.UserfunctionsInfo.label));
    gds.options.UserfunctionsInfo(dim+1).label = label;
    gds.options.UserfunctionsInfo(dim+1).name = name;
    gds.options.UserfunctionsInfo(dim+1).state = 1;
    gds.userfunction{dim+1} = get(handles.edituserfunction,'String');
    i = find(strcmp(namestr,name));
    if i
        button = questdlg('Userfunction already exist! Do you want to continue? If you press yes to continue, you will overwrite the existing userfunction',...
            'Userfunction already exist','Yes','No','No');
        if strcmp(button,'No')
            gds.userfunction(dim+1) = [];
            gds.options.UserfunctionsInfo(dim+1) = [];
            return;
        elseif strcmp(button,'Yes')
            gds.userfunction(dim+1)=[];
            gds.options.UserfunctionsInfo(dim+1) = [];
            gds.options.UserfunctionsInfo(i).label = label;
            gds.options.UserfunctionsInfo(i).name = name;
            gds.userfunction{i} = get(handles.edituserfunction,'String');
        end
    end
    if find(strcmp(labelstr,label))
        warndlg('There is another userfunction using this label: please choose another label!');
        set(handles.label,'String','');
        gds.userfunction(dim+1) = [];
        gds.options.UserfunctionsInfo(dim+1) = [];
    end
else
    gds.options.UserfunctionsInfo=[];
    gds.userfunction{1}=get(handles.edituserfunction,'String');
    gds.options.UserfunctionsInfo.label = label;
    gds.options.UserfunctionsInfo.name  = name;
    gds.options.UserfunctionsInfo.state = 1;
end
load_listbox1(handles);

  
% --------------------------------------------------------------------
function listbox1_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.renamebutton.
global gds;
if isfield(gds,'userfunction')
    val=get(handles.listbox1,'Value');
    d=0;
    for k=1:size(gds.userfunction,2)
        if ~isempty(gds.userfunction{k})
            d=d+1;
            if d==val
                set(handles.name,'String',gds.options.UserfunctionsInfo(k).name);
                name_Callback(handles,[],handles);
                set(handles.label,'String',gds.options.UserfunctionsInfo(k).label);
                set(handles.edituserfunction,'String',gds.userfunction{k});
                break;
            end
        end
    end
else
    set(handles.name,'String','');
    name_Callback(handles,[],handles);
    set(handles.label,'String','');
    set(handles.edituserfunction,'String','');
end


% --------------------------------------------------------------------
function label_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.renamebutton.
global MC
label=get(handles.label,'String');
if length(label)>2
    warndlg('A label exist of 2 characters')
    set(handles.label,'String',label(1:2));
elseif length(label)<2
    while length(label)<2
        label=sprintf('%s%c',label,char(32));
    end
    set(handles.label,'String',label(1:2));
end
list  = get(MC.mainwindow.initial_point,'children');
tag   = get(list,'Tag');
label = strcat(deblank(label),'_');
i = strmatch(label,tag);
if find(i)
    warndlg('Labels should be unique and differ from labels of standard special points.')
    set(handles.label,'String','');
end

% --------------------------------------------------------------------
function name_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.renamebutton.
global gds
set(handles.adbutton,'Enable','on');
set(handles.updatebutton,'Enable','off');
name = get(handles.name,'String');
if isfield(gds,'userfunction')%userfunctions already exists
    if find(strcmp(cellstr(strvcat(gds.options.UserfunctionsInfo.name)),name))
        set(handles.updatebutton,'Enable','on');
        set(handles.adbutton,'Enable','off');
    end
end


% --------------------------------------------------------------------
function deletebutton_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.deletebuttton.
global gds
% set(handles.deletebutton,'Enable','off');

name  = get(handles.listbox1,'String');
val   = get(handles.listbox1,'Value');
label = get(handles.label,'String');
if isfield(gds,'userfunction')%userfunctions already exists
    dim = size(gds.userfunction,2);
    i = find(strcmp(cellstr(char(gds.options.UserfunctionsInfo.name)),name{val}));
    if i
        button = questdlg('Do you want to continue? If you press yes to continue, you will delete the existing userfunction',...
                'Delete userfunction ','Yes','No','No');          
        if strcmp(button,'Yes')
            if isfield(gds,'poincare_eq') && ~isempty(gds.poincare_eq)
                funct=func2str(gds.poincare_eq);
                if strcmp(funct,name{val})
                    gds.poincare_eq=[];
                    gds.poincare_do=0;
                    gds=rmfield(gds,'poincare_eq');
                end
            end    
            gds.userfunction{i} = '';
            gds.options.UserfunctionsInfo(i).state = 0;
        end
    end
end
load_listbox1(handles);

% --------------------------------------------------------------------
function updatebutton_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.updatebutton.
global gds;
nameu = get(handles.name,'String'); label = get(handles.label,'String');
if isfield(gds,'userfunction')%userfunctions already exists
    i = find(strcmp(cellstr(char(gds.options.UserfunctionsInfo.name)),nameu));
    gds.options.UserfunctionsInfo(i).label = label;
    gds.options.UserfunctionsInfo(i).name = nameu;
    gds.options.UserfunctionsInfo(i).state = 1;
    gds.userfunction{i}=get(handles.edituserfunction,'String');
end
load_listbox1(handles);

% --------------------------------------------------------------------
function cancelbutton_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.cancelbutton.
global gds oldgds path_sys;
gds  = oldgds;
file = fullfile(path_sys,gds.system);
save(file,'gds');
delete(handles.userfunfig);
if gds.open.figuur==1,starter;end
if gds.open.numeric_fig==1,numeric;end



%-----------------------------------------------------------------
function load_listbox1(handles)
global gds;
if isfield(gds,'userfunction') && ~isempty(char(gds.userfunction))
    str = [];d=0;
    for i=1:size(gds.userfunction,2)
        if ~isempty(gds.userfunction{i})
            d=d+1;
            str{d,1} = gds.options.UserfunctionsInfo(i).name;
            val=i;
        end
    end
    set(handles.listbox1,'String',cellstr(str),'Value',d);
    set(handles.name,'String',gds.options.UserfunctionsInfo(val).name);
    name_Callback(handles,[],handles);
    set(handles.label,'String',gds.options.UserfunctionsInfo(val).label);
    set(handles.edituserfunction,'String',gds.userfunction{val});
else
    set(handles.listbox1,'String','');
    set(handles.name,'String','');
    name_Callback(handles,[],handles);
    set(handles.label,'String','');
    set(handles.edituserfunction,'String','res=');

end


% --------------------------------------------------------------------
function okbutton_Callback(h, eventdata, handles, varargin)
% Stub for Callback of the uicontrol handles.okbutton.
global gds path_sys MC;
name = get(handles.name,'String');
ad   = get(handles.adbutton,'Enable');
if ~isempty(name) && strcmp(ad,'on')
    adbutton_Callback(h,[],handles,[]);
end
updat = get(handles.updatebutton,'Enable');
if ~isempty(name) && strcmp (updat,'on')
    updatebutton_Callback(h,[],handles,[]);
end
load_matcont;
if isfield(gds,'userfunction'),gds.options=contset(gds.options,'Userfunctions',1);
else gds.options=contset(gds.options,'Userfunctions',0);end
string_jac = '';string_jacp = '';string_hess = '';string_hessp = '';string_tensor3='';string_tensor4='';string_tensor5='';
dimp = size(gds.parameters,1);par='';pa='';
if ~isempty(dimp)
    par = cellstr((strcat(',',char(gds.parameters{:,1}))));
    par = strcat(par{:,1});
    pa  = par(2:end);
end
cor = '';
t = gds.time{1,1};
if (~isempty(gds.dim))
    cor = cellstr((strcat(',',char(gds.coordinates{:,1}))));
    cor = strcat(cor{:,1}); cor = cor(2:end);
end
if (~isempty(t))
    t=strcat(t,',');
else
    t='t,';
end
string_sys = cellstr(systems('replace_sys_input',gds.equations));
if strcmp(string_sys,'error')
    errordlg('The number of equations is not correct','Error');
    return
end
fwrite = strcat(gds.system,'.m');
fwrite = fullfile(path_sys,fwrite);
[fid_write,message] = fopen(fwrite,'w');
if fid_write==-1
    errordlg(message,'Error');
    return
end
fread = fullfile(path_sys,'standard.m');
[fid_read,message] = fopen(fread,'r');
if fid_read == -1
    errordlg(message,'Error');
    return
end
string_handles={'out{1} = @init;';
                'out{2} = @fun_eval;';
                'out{3} = [];';
                'out{4} = [];';
                'out{5} = [];';
                'out{6} = [];';
                'out{7} = [];';
                'out{8} = [];';
                'out{9} = [];';
                'return;';};
string_init = cellstr(systems('make_init'));
if (gds.der(3,1)==1|gds.der(4,1)==1)
    str_init=string_init{2};
    str_handles1=string_handles{3};
    str_handles2=string_handles{4};
    str_on='''Jacobian'',handles(3)';str_off='''Jacobian'',[]';
    strp_on='''JacobianP'',handles(4)';strp_off='''JacobianP'',[]';
    if ~isempty(par)
        str_init=strrep(str_init,str_off,str_on);
    end
    str_handles1=strrep(str_handles1,'[]','@jacobian');
    string_handles{3,1}=str_handles1;
    str_handles2=strrep(str_handles2,'[]','@jacobianp');
    string_handles{4,1}=str_handles2;
    str_init=strrep(str_init,strp_off,strp_on);string_init{2,1}=str_init;
    string_jac=gds.jac;string_jacp=gds.jacp;
end
if (exist('sym')==2&& gds.der(4,1)==1)
    string_jac  = cellstr(systems('symjac',handles,gds.equations,cor,pa,10));
    string_jacp = cellstr(systems('symjac',handles,gds.equations,cor,pa,11));
end    
if gds.der(3,1)==1
    string_jac  = cellstr(systems('replace_jac_input',gds.jac));
    string_jacp = cellstr(systems('replace_jacp_input',gds.jacp));
end
if (gds.der(3,2)==1|gds.der(4,2)==1)
    str_init=string_init{2};
    str_handles1=string_handles{5};
    str_handles2=string_handles{6};
    str_on='''Hessians'',handles(5)';str_off='''Hessians'',[]';
    strp_on='''HessiansP'',handles(6)';strp_off='''HessiansP'',[]';
    if ~isempty(par)
    str_init=strrep(str_init,str_off,str_on);
    end
    str_handles1=strrep(str_handles1,'[]','@hessians');
    string_handles{5,1}=str_handles1;
    str_handles2=strrep(str_handles2,'[]','@hessiansp');
    string_handles{6,1}=str_handles2;
    str_init=strrep(str_init,strp_off,strp_on);string_init{2,1}=str_init;
    string_hess=gds.hess;string_hessp=gds.hessp;
end
if (exist('sym')==2&& gds.der(4,2)==1)
    string_hess  = cellstr(systems('symjac',handles,gds.equations,cor,pa,20));
    string_hessp = cellstr(systems('symjac',handles,gds.equations,cor,pa,21));
end    
if (gds.der(3,2)==1)
    string_hess  = cellstr(systems('replace_hess_input',gds.hess));
    string_hessp = cellstr(systems('replace_hessp_input',gds.hessp));
end
if (gds.der(4,3)==1), string_tensor3 = gds.tensor3;end
if (exist('sym')==2&& gds.der(4,3)==1)
    string_tensor3 = cellstr(systems('symjac',handles,gds.equations,cor,pa,3));
    str_handles2=string_handles{7};
    str_handles2=strrep(str_handles2,'[]','@der3');
    string_handles{7,1}=str_handles2;
end    
if (gds.der(4,4)==1), string_tensor4 = gds.tensor4;end
if (exist('sym')==2&& gds.der(4,4)==1)
    string_tensor4 = cellstr(systems('symjac',handles,gds.equations,cor,pa,4));
    str_handles2=string_handles{8};
    str_handles2=strrep(str_handles2,'[]','@der4');
    string_handles{8,1}=str_handles2;
end    
if (gds.der(4,5)==1), string_tensor5 = gds.tensor5;end
if (exist('sym')==2&& gds.der(4,5)==1)
    string_tensor5 = cellstr(systems('symjac',handles,gds.equations,cor,pa,5));
    str_handles2=string_handles{8};
    str_handles2=strrep(str_handles2,'[]','@der5');
    string_handles{8,1}=str_handles2;
end 
if ~isempty(gds.options.UserfunctionsInfo)
    siz = size(gds.options.UserfunctionsInfo,2);
    for i = 1:siz
        string_handles{9+i,1}= sprintf('out{%d}= @%s;',9+i,gds.options.UserfunctionsInfo(i).name);
    end
else siz=0;end
h=0;
while feof(fid_read)==0
    tline   = fgetl(fid_read);
    h=h+1;
    if h==2
        for i=1:9+siz
            fprintf(fid_write,'%s\n',string_handles{i,1});
        end
    end        
    matches=strrep(tline,'time,',t);
    matches=strrep(matches,'odefile',gds.system);
    matches=strrep(matches,',parameters',par);
    fprintf(fid_write,'%s\n',matches);
    if isfield(gds,'userfunction')
        if ~isempty(findstr(matches,'varargout{1}=der5(coordinates,'))          
            for i = 1:size(gds.userfunction,2)
                hs1 = sprintf('case ''%s''\n\tvarargout{1}=%s(coordinates%s);',gds.options.UserfunctionsInfo(i).name,gds.options.UserfunctionsInfo(i).name,par);
                fprintf(fid_write,'%s\n',hs1);
            end
        end
    end        
    if ~isempty(findstr(matches,'function dydt'))
        dim = size(string_sys,1);         
        for i=1:dim
              string_sys{i};
              fprintf(fid_write,'%s\n',string_sys{i});
        end
    end
    if ~isempty(findstr(matches,'function [tspan,y0,options]'))
        dim = size(string_init,1);  
        for i=1:dim
              fprintf(fid_write,'%s\n',string_init{i});
        end
    end
    if (~isempty(findstr(matches,'function jac '))&& ~isempty(string_jac))
        dim = size(string_jac,1);         
        for i = 1:dim
              fprintf(fid_write,'%s\n',string_jac{i});
        end
    end   
    if (~isempty(findstr(matches,'function jacp'))&& ~isempty(string_jacp))
       dim = size(string_jacp,1);       
       for i = 1:dim
           fprintf(fid_write,'%s\n',string_jacp{i});
       end
    end    
    if (~isempty(findstr(matches,'function hess '))&& ~isempty(string_hess))
        dim = size(string_hess,1);        
        for i=1:dim
              fprintf(fid_write,'%s\n',string_hess{i});
        end
    end
    if (~isempty(findstr(matches,'function hessp'))&& ~isempty(string_hessp))
        dim = size(string_hessp,1);
        for i=1:dim
              fprintf(fid_write,'%s\n',string_hessp{i});
        end
    end
    if (~isempty(findstr(matches,'function tens3'))&& ~isempty(string_tensor3))
        dim = size(string_tensor3,1);
        for i = 1:dim
              fprintf(fid_write,'%s\n',string_tensor3{i});
        end
    end
    if (~isempty(findstr(matches,'function tens4'))&& ~isempty(string_tensor4))
        dim = size(string_tensor4,1);
        for i=1:dim
              fprintf(fid_write,'%s\n',string_tensor4{i});
        end
    end
    if (~isempty(findstr(matches,'function tens5'))&& ~isempty(string_tensor5))
        dim = size(string_tensor5,1);
        for i = 1:dim
              fprintf(fid_write,'%s\n',string_tensor5{i});
        end
    end
end
if ~isempty(gds.options.UserfunctionsInfo)    
   for i=1:size(gds.options.UserfunctionsInfo,2)
       res=0;
       if isfield(gds,'userfunction') && ~isempty(gds.userfunction{i})
           str_user = systems('replace_token',cellstr(gds.userfunction{i}));
       else 
           str_user=cellstr('res=');
       end
       hs1 = sprintf('function userfun%d=%s(t,kmrgd%s)',i,gds.options.UserfunctionsInfo(i).name,par);
       fprintf(fid_write,'%s\n',hs1);
       hs1 = sprintf('userfun%d',i);
       dim = size(str_user,1);
       for j = 1:dim
           d = strmatch('res=',str_user{j},'exact');
           if findstr('res',str_user{j}),res=1;end
           str_user{j} = strrep(str_user{j},'res',hs1);
           if d==1
               fprintf(fid_write,'\t%s=0;\n',hs1);
           else 
               fprintf(fid_write,'\t%s;\n',str_user{j});
           end
       end
       if res==0,fprintf(fid_write,'\t%s=0;\n',hs1);end
   end
end        
fclose(fid_read);
fclose(fid_write);
file = fullfile(path_sys,gds.system);
delete(handles.userfunfig);
save(file,'gds');
if gds.open.figuur==1,starter;end
if gds.open.numeric_fig==1,numeric;end
rehash;

