function gui_DO
global gds
    gds.gui_point = @point;
    gds.gui_type = @curvetype;
    gds.gui_starter = @starter1;
    gds.gui_load_layout = @load_layout;
    gds.gui_layoutbox = @layoutbox;
    gds.gui_numeric = @numeric1;
    gds.gui_load_draw = @load_draw;
    gds.gui_load_attributes = @load_attributes;
    gds.gui_make_ready = @make_ready;
    gds.gui_label = @label1;
    gds.gui_numeric_label = @numeric_label;
    gds.gui_draw = @draw;
    gds.gui_output = @output;
    gds.gui_load_point = @load_point;
    gds.gui_singularities = @singularities;
    gds.gui_start_cont = @start_cont;

%-------------------------------------------------------------------------
function point(tag)        
global gds
    [point,str]=strtok(tag,'_');
    set(0,'ShowHiddenHandles','on');
    h=findobj('Type','uimenu','Tag','window');
    set(h,'enable','on');
    h=findobj('Type','uimenu','Tag','curvetype');
    set(0,'ShowHiddenHandles','off');
    list = get(h,'children');
    type = get(list,'tag');
    str  = strcat(str,'_');
    for i=1:length(list)
        if isempty(findstr(str,type{i}))
            set(list(i),'enable','off');
        else
            set(list(i),'enable','on');
        end
    end     
    gds.point=sprintf('%s%c',point,char(32));         
   curvetype;%callback of 'Type->curve->equilibrium'
%-------------------------------------------------------------------------    
function curvetype
global gds path_sys MC
    if isempty(gds.point)
        gds.point='P ';
    end
    gds.type='DO ';
    matcont('make_curve');
    load_matcont;
    file=strcat(gds.system,'.mat');
    file=fullfile(path_sys,file);
    save(file,'gds');
    integrator;starter;
    if  ~isempty(MC.numeric_fig),numeric;end
    set(0,'ShowHiddenHandles','on');        
    h=findobj('Type','uimenu','Tag','window');
    set(h,'enable','on');
    set(0,'ShowHiddenHandles','off');     
%-------------------------------------------------------------------------   
function starter1(handles)
global gds 
    ndim=size(gds.parameters,1);
    color = [1 1 1];
    gds.options.ActiveParams=[];
    gds.options.ActiveUParams=[];
    gds.options.ActiveSParams=[];
    gds.options.ActiveSParam=[];    
    s = gds.dim*2+ndim*2+12;
    slider_step(1) = 2/s;
    slider_step(2) = 2/s;
    set(handles.slider1,'sliderstep',slider_step,'max',s,'min',0,'Value',s);
    j = s-2;
    stat=uicontrol(handles.figuur,'Style','text','String','Initial Point','Tag','Initial_Point','BackGroundColor',color,'Units','characters','fontname','FixedWidth','fontsize',12);
    pos = [10 j 30 1.8];user.num = 0;user.pos = pos;
    set(stat,'Position',pos,'UserData',user);
    pos1 = starter('start_time',handles,j);
    pos1 = starter('start_coordinates',handles,pos1);
    pos1 = starter('start_parameters_orbit',handles,pos1); 
    start_PMfunction(handles,pos1);
    starter('in_start',handles);
%-------------------------------------------------------------------------   
function  string = load_layout(handles)
global gds
    for i=1:3
        if gds.numeric.O{i,2}==1
            gds.numeric.O{i,1}=upper(gds.numeric.O{i,1});
        elseif gds.numeric.O{i,2}==0
            gds.numeric.O{i,1}=lower(gds.numeric.O{i,1});
        end
        string{i,1}=gds.numeric.O{i,1};
    end    
    set(handles.layoutbox,'String',string);
%-------------------------------------------------------------------------
function layoutbox(list,index_selected)
global gds
    d=gds.numeric.O{index_selected,2};    
    gds.numeric.O{index_selected,2}=1-d;    

%-------------------------------------------------------------------------    
function numeric1(handles)
global gds 
    ndim=size(gds.parameters,1);
    s = 2;
    if gds.numeric.O{1,2}==1
        s = s+3;
    end
    if gds.numeric.O{2,2}==1
        s = s+2*gds.dim+2;
    end 
    if gds.numeric.O{3,2}==1        
        s = s+2*ndim+2;
    end     
    slider_step(1) = 2/s;
    slider_step(2) = 2/s;
    set(handles.slider1,'sliderstep',slider_step,'max',s,'min',0,'Value',s);
    j=s;    
    if gds.numeric.O{1,2}==1
        j=numeric('start_time',handles,j);
    end
    if gds.numeric.O{2,2}==1
        j=numeric('start_coordinates',handles,j);
    end 
    if gds.numeric.O{3,2}==1        
        j=numeric('start_parameters',handles,j);
    end     
    numeric('in_numeric',handles);
 %-------------------------------------------------------------------------   
function [string,plotopts] = load_draw
global gds
    plotopts = {};
    string={'Coordinates';'Parameters';'Time'};
    if ~isempty(gds.options.Userfunctions) && gds.options.Userfunctions ~= 0
        string{end+1} = 'Userfunction';
    end
%-------------------------------------------------------------------------
function e = make_ready(plo,x)
    len=size(plo,2);
    e(1,len) = 0;
    for j=1:len
        switch plo(j).type
        case 'Coordinates'
            e(1,j)=plo(j).val;
        case 'Parameters'
            e(1,j)=-plo(j).val;
        case 'Time'
            e(1,j)=0;
        case 'Userfunction'
            e(1,j) = -100000*ndim - pl(j).val;
        otherwise
            for j = 1:len
                e(1,j)=inf;
                return
            end
        end 
    end    
%------------------------------------------------------
function label = numeric_label(numsing) 
global gds MC
    num = 1;
    label{gds.dim} = '';
    if isfield(MC.numeric_handles,'coord1')
        for d=1:gds.dim
            label{num}=sprintf('set(MC.numeric_handles.coord%d,''String'',num2str(x(%d,i),''%%.8g''))',d,d);
            num = num+1;
        end
    end
    if isfield(MC.numeric_handles,'time')
        label{num}=sprintf('set(MC.numeric_handles.time,''String'',num2str(t(i,1),''%%.10g''))'); 
        num = num+1;
    end
    if isfield(MC.numeric_handles,'param1')
        for d=1:size(gds.parameters,1);
            label{num}=sprintf('set(MC.numeric_handles.param%d,''String'',gds.parameters{%d,2})',d,d);
            num = num+1;
        end 
    end
%-------------------------------------------------------------------------    
function draw(varargin)
global gds MC
    ndim = size(gds.parameters,1);
    if strcmp(varargin{4},'select')      
        if nargin < 7
            x = varargin{3}; s = varargin{5};t = varargin{2}; 
        else
            x = varargin{3}; s = varargin{5};h = varargin{6}; f = varargin{7};
        end
    else
        file = varargin{3};load(file);
        s(1)   = [];   s(end) = [];
    end
    switch varargin{1}
    case 2
        plot2=varargin{2};
        if (plot2==[inf inf]),return;end
        hold on;d=axis;
        p=0;axis(d);
        plo2{1}='empt';plo2{2}='empt' ;
        k=size(x,2);
        skew = 0.01*[d(2)-d(1) d(4)-d(3)];
        axis(d);plo2{1} = 'empt';plo2{2} = 'empt' ;    
        plo2s{1} = 'empt';plo2s{2} = 'empt';      
        p = size(cat(1,s.index),1);
        if p>0,        sdata = cat(1,s.data); end
        for j=1:2
            if plot2(1,j) <= -100000*ndim % case Userfunction
                plo2{j}  = h(2-plot2(1,j)-100000*ndim,:);
            
            elseif plot2(1,j)<0
                pl2=-plot2(1,j);
                plo2{j}=gds.parameters{pl2,2}*ones(1,k);
                if p>0, plo2s{j} = gds.parameters{pl2,2}*ones(1,p);end
            end      
            if plot2(1,j)==0
                plo2{j}=t(1:k);
                if p>0, plo2s{j}=cat(1,sdata.t);end
            end
            if strcmp(plo2{j},'empt')==1
                plo2{j}=x(plot2(1,j),1:k);
                if p>0
                    sdatax=cat(1,sdata.x)';
                    plo2s{j} =  sdatax(plot2(1,j),:);
                end
            end
        end  
        if strcmp(varargin{4},'redraw')      
            if ~strcmp(plo2{1},'empt')&&~strcmp(plo2{2},'empt')
              line(plo2{1},plo2{2},'LineStyle','-','Color','b');
            end
           if ~strcmp(plo2s{1},'empt')&~strcmp(plo2s{2},'empt') & p>0
               line(plo2s{1},plo2s{2},'LineStyle','none', 'Marker','*','Color','r');
           end     
        else
            plot(plo2s{1},plo2s{2},'m+','MarkerSize',200);                
             hold off;
        end
    case 3 %3D
        plo=varargin{2};
        if (plo==[inf inf inf]),return;end
        hold on;d=axis;
        p=0;axis(d);
        skew = 0.01*[d(2)-d(1) d(4)-d(3) d(6)-d(5)];
        k = size(x,2);p=size(cat(1,s.index),1);
         if p>0,        sdata = cat(1,s.data);end
        axis(d);plo3{1} = 'empt';plo3{2} = 'empt';plo3{3} = 'empt';
        plo3s{1} = 'empt'; plo3s{2} = 'empt';plo3s{3} = 'empt';
        k=size(x,2);
        for j=1:3
            if plo(1,j) <= -100000*ndim % case Userfunction
                plo3{j}  = h(2-plo(1,j)-100000*ndim,:);
            elseif plo(1,j)<0
                pl3=-plo(1,j);
                plo3{j}  = gds.parameters{pl3,2}*ones(1,k);
                if p>0, plo3s{j} = gds.parameters{pl2,2}*ones(1,p);end
            end               
            if plo(1,j)==0
                plo3{j}  = t(1:k);
                if p>0,plo3s{j} = cat(1,sdata.t);end
            end
            if strcmp(plo3{j},'empt')==1
                plo3{j}  = x(plo(1,j),1:k);
                if p>0
                    sdatax   = cat(1,sdata.x)';
                    plo3s{j} =  sdatax(plo(1,j),:);
                end
            end
        end  
        if strcmp(varargin{4},'redraw')      
            if ~strcmp(plo3{1},'empt')&&~strcmp(plo3{2},'empt')&&~strcmp(plo3{3},'empt')
                line(plo3{1},plo3{2},plo3{3},'linestyle','-','Color','b');
            end
             if ~strcmp(plo3s{1},'empt')&~strcmp(plo3s{2},'empt')&~strcmp(plo3s{3},'empt')& p>0
                line(plo3s{1},plo3s{2},plo3s{3},'linestyle','none' , 'Marker', '*','Color','r');
             end
        else
            plot3(plo3s{1},plo3s{2},plo3s{3},'m+','MarkerSize',200);                
            hold off;
        end
    case 4 %numeric
        for d=1:ndim
            parameters(d) = gds.parameters{d,2};
        end   
        dat = get(MC.numeric_fig,'Userdata');
        i = s.index;
        for k=1:size(dat.label,2)
            eval(dat.label{k});  
        end 
    end  
%-------------------------------------------------------------------------    
function output
%-------------------------------------------------------------------------   
function load_point(index,x,string,file,num)
global gds
    load(file);
    if isfield(s(num).data,'x')
        gds.time{1,2}=s(num).data.t;
        for i=1:gds.dim
            gds.coordinates{i,2}=s(num).data.x(i);
        end
    else
        gds.time{1,2}=t(index);
        for i=1:gds.dim
            gds.coordinates{i,2}=x(i,index);
        end
    end
    val=gds.integrator.tspan(2)-gds.integrator.tspan(1);
    gds.integrator.tspan=[gds.time{1,2} (val+gds.time{1,2})];
    if exist('option','var')
        gds.integrator.options = option;
    else
        option = [];
    end
    ndim=size(param,1);
    for i=1:ndim
        gds.parameters{i,2}=param(i,1);
    end
    if strcmp(string,'save')
        save(file,'x','s','param','t','ctype','point','option');
    end
%-------------------------------------------------------------------------    
function start_cont(varargin)
global gds path_sys MC sys sp
    ndim=size(gds.parameters,1);
    if ~isempty(varargin)
        funhandle = feval(gds.system);
        file = fullfile(path_sys,gds.system,gds.diagram,strcat(gds.curve.new,'.mat'));
        if exist(file,'file')
            load(file);
        else
            errordlg('It is not possible to extend the current curve!','Error extend EP');
            return 
        end  
        x0 = x(:,end);
        for k = 1:ndim
            gds.parameters{k,2} = param(k);
        end
        if exist('option','var')
            gds.integrator.options = option;
        end
        set(0,'ShowHiddenHandles','on');
        hh = findobj('Style','edit','tag','interval');
        val = str2double(get(hh,'String'));
        gds.time{1,2} = abs(t(end));
        if val<0
            errordlg('this isn''t possible');
            set(hh,'string','1');
        else
            gds.integrator.tspan=[gds.time{1,2} (val+gds.time{1,2})];
        end
    else %case extend 
        funhandle = feval(gds.system);
        x0 = vertcat(gds.coordinates{:,2});
        gds.integrator.options = odeset(gds.integrator.options,'Jacobian',funhandle{3});
        gds.integrator.options = odeset(gds.integrator.options,'Hessians',funhandle{5});
    end
    try
        if (gds.options.Backward==1)
            interval = abs(gds.integrator.tspan(2)-gds.integrator.tspan(1));
            tspan = [gds.integrator.tspan(1) gds.integrator.tspan(1)-interval];
        else
            tspan = gds.integrator.tspan;
        end
        gds.integrator_fun=funhandle{2};
        str = sprintf('(funhandle{2},tspan,x0,gds.integrator.options');
        par = cellstr(strcat(',',num2str(vertcat(gds.parameters{:,2}))));
        gds.integrator.options = odeset(gds.integrator.options,'OutputFcn',@integ_plotDO);
        str = strcat(str,strcat(par{:,1}),')');
        switch gds.integrator.method
        case 'ode45'
            str = strcat('ode45',str);
        case 'ode23'
            str = strcat('ode23',str);
        case 'ode113'
            str = strcat('ode113',str);
        case 'ode15s'
            str = strcat('ode15s',str);
        case 'ode23s'
            str = strcat('ode23s',str);
        case 'ode23t'
            str = strcat('ode23t',str);
        case 'ode23tb'
            str = strcat('ode23tb',str);
        end  
        npoints=1;
        npoints = sys.gui.plot_points;
        gds.integrator.options = odeset(gds.integrator.options,'OutputFcn',@integ_plotDO);
        file = fullfile(path_sys,gds.system);
        save(file,'gds');
        set(0,'ShowHiddenHandles','on');
        duration = findobj('style','text','Tag','duration');
        status = findobj('style','text','tag','mstatus');
        set(0,'ShowHiddenHandles','off');
        set(status,'string','computing');
        StartTime = clock;
        sp(1,1).index = 1;         sp(1,1).msg = 'This is the first point of the orbit'; sp(1,1).label = '00';
        sp(1,1).data=[];
        h=0;dif=[];k=0;
        [t,y] = eval(str);
        EndTime = clock;
        string  = sprintf('%.1f secs\n', etime(EndTime, StartTime));
        set(duration,'String',string);
        set(status,'String','ready');
        x = y';ind1=size(x,2);
        file = fullfile(path_sys,gds.system,gds.diagram,gds.curve.new);
        if isempty(varargin)
            sp(1,1).data.t=t(1);sp(1,1).data.x=x(:,1)';
            sp(end+1,1).index =ind1; sp(end,1).msg = 'This is the last point of the orbit';  sp(end,1).label = '99';sp(end,1).data=[];
            sp(end,1).data.t=t(end);sp(end,1).data.x=x(:,end)';
            s=sp;            
        else
            x1 = x;t1 = t;
            load(file);
             s(end)=[];%remove last point on the original curve
            x = [x,x1];
            t = [t;t1];
            ind1=size(x,2);
            sp(1)=[];
            s=[s;sp];s(1,1).data.t=t(1);s(1,1).data.x=x(:,1)';
            s(end+1,1).index = ind1; s(end,1).msg = 'This is the last point of the orbit';  s(end,1).label = '99';s(end,1).data=[];
            s(end,1).data.t=t(end);s(end,1).data.x=x(:,end)';
        end  
        ctype  = gds.type; point = gds.point;        
        status = mkdir(path_sys,gds.system);
        dir = fullfile(path_sys,gds.system);
        status = mkdir(dir,gds.diagram);
        param  = vertcat(gds.parameters{:,2});  
        option = gds.integrator.options;
        save(file,'t','x','s','param','ctype','point','option'); 
        sp=[];
        clear sp;
        file=fullfile(path_sys,gds.system);save(file,'gds');
        sys.gui.plot_points = npoints;
    catch
        errordlg(lasterr,'error integration');
    end 


%---------------------------------------------------------------------------------------
function j = start_select_cycle(handles,j)
%enter field for select cycle button
global gds path_sys;
j = j-3;
pos  = [10 j 30 2]; user.num = 0; user.pos = pos;
stat = uicontrol(handles.figuur,'Style','pushbutton','String','Select Cycle','Tag','select_cycle','units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',pos,'UserData',user,'Callback','selectcycle');
guidata(handles.figuur,handles);
%---------------------------------------------------------------------------------------
function j = start_PMfunction(handles,j)
%enter field for amplitude
global gds path_sys;
if ~isfield(gds,'poincare_do')||~isfield(gds,'poincare_eq'), gds.poincare_do=1;gds.poincare_eq='';gds.integrator.Crossection='';gds.integrator.nr_crossection=0;end
color = [1 1 1];
j = j-2;
user.num = 0; 
post = [2 j 18 1.80];user.pos=post;
pose = [20 j 25 1.80];
stat = uicontrol(handles.figuur,'Style','text','HorizontalAlignment','left','String','Crosssection','Tag','text_DOfunction','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',post,'UserData',user);
index=1; string{1}='none'; string{2}='add function';
if isfield(gds,'userfunction')&& ~isempty(char(gds.userfunction))
    dimu = size(gds.userfunction,2);
    d=0;
    for k= 1:dimu
        if ~isempty(gds.userfunction{k})
            d=d+1;
            string{2+d} = gds.options.UserfunctionsInfo(k).name;
            if strcmp(gds.options.UserfunctionsInfo(k).name,gds.integrator.Crossection)
                index=2+d;
            end
        end
    end
end
edit = uicontrol(handles.figuur,'Style','popupmenu','HorizontalAlignment','left','String',string,'Value',index,'Tag','DOfunction','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
try gds.poincare_eq=str2func(gds.poincare_eq);catch end
user.pos = pose; set(edit,'Position',pose,'Callback','set_option','UserData',user);
j = j-2;
post = [2 j 18 1.80];user.pos=post;
pose = [20 j 25 1.80];
stat = uicontrol(handles.figuur,'Style','text','HorizontalAlignment','left','String','Direction','Tag','nr_DOfunction','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',post,'UserData',user);
edit = uicontrol(handles.figuur,'Style','popupmenu','String',{'all','increasing','decreasing'},'HorizontalAlignment','left','Tag','nr_DOfunction','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
user.pos = pose; set(edit,'Position',pose,'Callback','set_option','UserData',user);

guidata(handles.figuur,handles);


