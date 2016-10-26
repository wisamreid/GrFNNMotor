function gui_HTHSN
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
global gds MC
[point,str]=strtok(tag,'_');
set(0,'ShowHiddenHandles','on');
h = findobj('Type','uimenu','Tag','window');
set(h,'enable','on');
h = findobj('Type','uimenu','Tag','curvetype');
set(0,'ShowHiddenHandles','off');
list = get(h,'children');
type = get(list,'tag');
str  = strcat(str,'_');
for i=1:length(list)
    if isempty(findstr(str,strcat(type{i},'_')))
        set(list(i),'enable','off');
    else
        set(list(i),'enable','on');
    end
end    
gds.point = sprintf('%s%c',point,char(32));
if strmatch('__',str,'exact')
    set(MC.mainwindow.compute,'enable','off');    
    gds.type='';   
    gds.curve.new='';
else
    curvetype;    
end

%-------------------------------------------------------------------------
function curvetype
global gds path_sys MC
if isempty(gds.point)
    gds.point = 'HTHSN ';
end
gds.type = 'HTHSN ';
set(0,'ShowHiddenHandles','on');        
h = findobj('Type','uimenu','Tag','forward'); set(h,'enable','on');
h = findobj('Type','uimenu','Tag','backward'); set(h,'enable','on');
h = findobj('Type','uimenu','Tag','extend'); set(h,'enable','on');
set(0,'ShowHiddenHandles','off');      
matcont('make_curve');
load_matcont;
file = fullfile(path_sys,strcat(gds.system,'.mat'));
save(file,'gds');
starter;
continuer; 
if  ~isempty(MC.numeric_fig),numeric;end
set(0,'ShowHiddenHandles','on');        
h = findobj('Type','uimenu','Tag','window');set(h,'enable','on');
h = findobj('Type','uimenu','Tag','extend'); set(h,'enable','off');
set(0,'ShowHiddenHandles','off');
    
%-------------------------------------------------------------------------
function starter1(handles)
global gds

ndim = size(gds.parameters,1);
color = [1 1 1];
gds.options.ActiveParams=[];
gds.options.ActiveUParams=[];
gds.options.ActiveSParams=[];
gds.options.ActiveT = [];
gds.options.ActiveEps1 = [];
gds.eps1tol = [];
gds.extravec = [0 0 0];
s = 2*gds.dim+ndim*2+42;
if isfield(gds,'userfunction')
    s = s+size(gds.userfunction,2)*2+2;   
end
slider_step(1) = 2/s;
slider_step(2) = 2/s;
set(handles.slider1,'sliderstep',slider_step,'max',s,'min',0,'Value',s);
j = s-2;
stat = uicontrol(handles.figuur,'Style','text','String','Initial Point','Tag','Initial_Point','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
pos = [10 j 35 1.80]; user.num = 0; user.pos = pos;
set(stat,'Position',pos,'UserData',user);

pos1 = starter('start_parameters',handles,j);
pos1 = start_USparams(handles,pos1);
pos1 = start_period(handles,pos1);
pos1 = start_eps1tol(handles,pos1);
%pos1 = start_eps1(handles,pos1);%veldje voor eps1
pos1 = starter('start_jacobian',handles,pos1);
pos1 = starter('start_discret',handles,pos1);
if isfield(gds,'userfunction')
    starter('start_userfunctions',handles,pos1);   
end           

starter('in_start',handles);
%-------------------------------------------------------------------------        

function load_layout(handles)
global gds 
  
if size(gds.numeric.HTHSN,1)<7
    gds.numeric.HTHSN{6,1} = 'user functions';    
    gds.numeric.HTHSN{6,2} = 0;   
    gds.numeric.HTHSN{7,1} = 'npoints';    
    gds.numeric.HTHSN{7,2} = 0;    
end
for i = 1:7
    if gds.numeric.HTHSN{i,2}==1    
        gds.numeric.HTHSN{i,1} = upper(gds.numeric.HTHSN{i,1});        
    elseif gds.numeric.HTHSN{i,2}==0   
        gds.numeric.HTHSN{i,1} = lower(gds.numeric.HTHSN{i,1});       
    end
    string{i,1} = gds.numeric.HTHSN{i,1};   
end
set(handles.layoutbox,'String',string);

%-------------------------------------------------------------------------

function layoutbox(list,index_selected)
global gds 

c = gds.numeric.HTHSN{index_selected,2};    
gds.numeric.HTHSN{index_selected,2} = 1-c; 

%-------------------------------------------------------------------------

function numeric1(handles)
global gds

ndim = size(gds.parameters,1);
if ~isfield(gds.numeric,'HTHSN') 
    gds.numeric.HTHSN = {'parameters' 1;'period' 0;'testfunctions' 1;'multipliers' 0;'current stepsize' 0};    
end

if size(gds.numeric.HTHSN,1)<6
    gds.numeric.HTHSN{6,1} = 'user functions';    
    gds.numeric.HTHSN{6,2} = 0;    
    gds.numeric.HTHSN{7,1} = 'npoints';    
    gds.numeric.HTHSN{7,2} = 0;    
end
gds.numeric.HTHSN{8,1} = 'eps1';
gds.numeric.HTHSN{8,2} = 1;    

for i = 1:size(gds.UParams,1)
    gds.numeric.HTHSN{8+i,1} = strcat('UParam',i);   
    gds.numeric.HTHSN{8+i,2} = 1;   
end
val = 8+size(gds.UParams,1);
for i = 1:size(gds.SParams,1)
    gds.numeric.HTHSN{val+i,1} = strcat('SParam',i);   
    gds.numeric.HTHSN{val+i,2} = 1;    
end

if isfield(gds,'userfunction')
    dimu = size(gds.userfunction,2);    
else dimu=0; 
end
s = 20;
dims = size(gds.SParams,1);
    
if gds.numeric.HTHSN{1,2}==1
    s = s+2*ndim+2;   
end
if gds.numeric.HTHSN{2,2}==1
    s = s+4;    
end
if gds.numeric.HTHSN{3,2}==1
    s = s+2*dims+2;   
end
if gds.numeric.HTHSN{4,2}==1
    s = s+4*gds.dim+2;   
end
if gds.numeric.HTHSN{5,2}==1
    s = s+2;    
end
if gds.numeric.HTHSN{6,2}==1
    s = s+2*dimu+2;   
end
if gds.numeric.HTHSN{7,2}==1
    s = s+2;   
end
if gds.numeric.HTHSN{8,2}==1
    s = s+2;    
end
for k = 1:size(gds.UParams,1)
    if gds.numeric.HTHSN{8+k,2}==1   
        s = s+2;   
    end
end
val = 8+size(gds.UParams,1);
for k = 1:size(gds.SParams,1)
    if gds.numeric.HTHSN{val+k,2}==1   
        s = s+2;        
    end
end

slider_step(1) = 2/s;
slider_step(2) = 2/s;
set(handles.slider1,'sliderstep',slider_step,'max',s,'min',0,'Value',s);
j = s;   
if gds.numeric.HTHSN{1,2}==1
    j = numeric('start_parameters',handles,j);    
end
if gds.numeric.HTHSN{3,2}==1
    j = start_testfunctions(handles,j);    
end
if gds.numeric.HTHSN{4,2}==1
    j = numeric('start_eigenvalues',handles,j);   
end
if gds.numeric.HTHSN{5,2}==1
    j = numeric('start_stepsize',handles,j);    
end
if gds.numeric.HTHSN{6,2}==1
    j = numeric('start_userfunctions',handles,j);    
end
if gds.numeric.HTHSN{7,2}==1
    j = numeric('start_npoints',handles,j);   
end
if gds.numeric.HTHSN{8,2}==1
    j = numeric_HTeps1(handles,j);   
end
if gds.numeric.HTHSN{9,2}==1
    numeric_HTpars(handles,j);        
end

numeric('in_numeric',handles);


%-------------------------------------------------------------------------    
function [string,plotopts] = load_draw
global gds
    
plotopts = {};
string = {'Coordinates';'Parameters';'Unstable Parameters';'Stable Parameters';'Period';'eps1';};
if gds.options.Eigenvalues==1
    string{7} = 'Eigenvalues';    
end
if ~isempty(gds.options.Userfunctions) && gds.options.Userfunctions ~= 0
    string{end+1} = 'Userfunction';    
end    
%-------------------------------------------------------------------------

function string = load_attributes
global gds 
  
d = 1;
for k = 1:gds.dim
    string{d,1} = sprinf('Mod[%d]',k);    
    string{d+1,1} = sprintf('Arg[%d]',k);   
    d = d+2;    
end
%-------------------------------------------------------------------------

function e = make_ready(pl,x)
global gds HTHSNds

ndim = size(gds.parameters,1);
len  = size(pl,2);
leng = size(x,1);
e(1,len) = 0;
    
for j = 1:len
    switch pl(j).type    
        case 'Coordinates'        
            e(1,j) = pl(j).val;                    
        case 'Parameters'
            e(1,j) = -pl(j).val;            
        case 'Unstable Parameters'           
            val = 0;            
            for k = 1:size(gds.options.ActiveUParams,2)            
                if pl(j).val==gds.options.ActiveUParams(1,k)
                    e(1,j) = HTHSNds.ncoords+k;
                    val = 1;
                end
            end
            if val == 0
                e(1,j) = -ndim-pl(j).val;
            end
            
        case 'Stable Parameters'              
            val = 0;            
            for k = 1:size(gds.options.ActiveSParams,2)            
                if pl(j).val==gds.options.ActiveSParams(1,k)                
                    e(1,j) = HTHSNds.ncoords+size(gds.options.ActiveUParams,2)+k;                    
                    val = 1;                    
                end
            end
            if val ==0            
                e(1,j) = -ndim-size(gds.UParams,1)-pl(j).val;                
            end
        case 'Period'
            e(1,j) = HTHSNds.PeriodIdx+1;
            
        case 'eps1'
            if HTHSNds.index == 1
                e(1,j) = HTHSNds.ncoords+size(gds.options.ActiveUParams,2)+size(gds.options.ActiveSParams,2)+1;
            else
                e(1,j) = HTHSNds.PeriodIdx+2;
            end
                
        case 'Testfunction'
            numsing = sing_numeric(gds.options.IgnoreSingularity);
            if (gds.options.Userfunctions==1)
                dimu = size(gds.userfunction,2);
            else dimu = 0;
            end
            e(1,j) = numsing(pl(j).val)+leng+dimu;    
            
        case 'Eigenvalues'
            if gds.options.Eigenvalues==0
                errordlg('It isn''t possible to plot the eigenvalues(they are not being computed)!');
            else
                e(1,j) = -ndim-size(gds.UParams,1)-size(gds.SParams,1)-pl(j).val;
            end
        
        case 'Userfunction'
            e(1,j) = -100000*ndim - pl(j).val;
        otherwise
            for j = 1:len            
                e(1,j)=inf;                
                return
            end
    end
end


%-------------------------------------------------------------------------

function [label,lab] = label1(plot2)
global gds

ndim = size(gds.parameters,1);
len = size(plot2,2); dimplot = 0;
label{len} = '';
lab{len} = '';
for k = 1:len
    label{k}='empty';    
    lab{k}='empty';    
    if plot2(1,k)<0   
        if (plot2(1,k)>=-ndim)
            pl2  = -plot2(1,k);        
            label{k}  = sprintf('gds.parameters{%d,2}*ones(HTHSNds.tps,p)',pl2);          
            lab{k}  = sprintf('gds.parameters{%d,2}',pl2);                      
        elseif (plot2(1,k)<-ndim) && (plot2(1,k)>=-ndim-size(gds.UParams,1))
            pl2  = -plot2(1,k)-ndim;            
            label{k}  = sprintf('gds.UParams{%d,2}*ones(HTHSNds.tps,p)',pl2);              
            lab{k}  = sprintf('gds.UParams{%d,2}',pl2);                  
        elseif (plot2(1,k)<-ndim-size(gds.UParams,1)) && (plot2(1,k)>=-ndim-size(gds.UParams,1)-size(gds.SParams,1))        
            pl2  = -plot2(1,k)-ndim-size(gds.UParams,1);           
            label{k}  = sprintf('gds.SParams{%d,2}*ones(HTHSNds.tps,p)',pl2);              
            lab{k}  = sprintf('gds.SParams{%d,2}',pl2);           
        end
    end    
    if plot2(1,k)==0    
        label{k}  = 'gds.time{1,2}*ones(HTHSNds.tps,p)';        
        lab{k} = 'gds.time{1,2}';        
    end

    if plot2(1,k) < -(ndim+size(gds.UParams,1)+size(gds.SParams,1)) && plot2(1,k) >-10000*ndim%case eigenvalues                   
        if mod(plot2(1,k)+ndim+size(gds.UParams,1)+size(gds.SParams,1),2)==0        
            label{k}  = sprintf('angle(fout(HTHSNds.ntst+1+%d,i))',(-plot2(1,k)-(ndim+size(gds.UParams,1)+size(gds.SParams,1)))/2);            
            lab{k} = sprintf('angle(fout(HTHSNds.ntst+1+%d,i))',(-plot2(1,k)-(ndim+size(gds.UParams,1)+size(gds.SParams,1)))/2);            
        else %Re        
            label{k} = sprintf('abs(fout(HTHSNds.ntst+1+%d,i))',(-plot2(1,k)-(ndim+size(gds.UParams,1)+size(gds.SParams,1))+1)/2);            
            lab{k} = sprintf('abs(fout(HTHSNds.ntst+1+%d,i))',(-plot2(1,k)-(ndim+size(gds.UParams,1)+size(gds.SParams,1))+1)/2);                    
        end
    elseif plot2(1,k) <= -100000*ndim % case Userfunction    
        label{k}  = sprintf('hout(2+%d,i)',-plot2(1,k)-100000*ndim);         
        lab{k}  = sprintf('hout(2+%d,i)',-plot2(1,k)-100000*ndim);         
    end
    
    if strcmp(label{k},'empty')==1    
        dimplot = dimplot+1;        
        if plot2(1,k) <= gds.dim %case coordinate        
            label{k}  = sprintf('xout((0:HTHSNds.tps-1)*HTHSNds.nphase+%d,i)',plot2(1,k));            
            lab{k}  = sprintf('xout((0:HTHSNds.tps-1)*HTHSNds.nphase+%d,i)]',plot2(1,k));                      
        else
            label{k} = sprintf('ones(HTHSNds.tps,1)*xout(%d,i)',plot2(1,k));            
            lab{k}  = sprintf('xout(%d,i)',plot2(1,k));             
        end
    end
end
if dimplot == 0
    for k = 1:len    
        label{k}=lab{k};        
    end
end
    
%-------------------------------------------------------------------------    
function label = numeric_label(numsing) 
global gds MC HTHSNds
 
ndim = size(gds.parameters,1);
num = 1;
label{ndim} = '';    
if isfield(MC.numeric_handles,'param1')
    for d = 1:ndim    
        label{num} = sprintf('set(MC.numeric_handles.param%d,''String'',num2str(parameters(%d),''%%0.8g''))',d,d);        
        num = num +1;        
    end
end

    if isfield(MC.numeric_handles,'UParam1')
        for k = 1:size(gds.UParams,1)    
            label{num} = sprintf('set(MC.numeric_handles.UParam%k,''String'',num2str(UParams(%k),''%%0.8g''))',k,k);        
            num = num+1;        
        end
    end
       
    if isfield(MC.numeric_handles,'SParam1')
            label{num} = sprintf('set(MC.numeric_handles.SParam1,''String'',num2str(SParams(1),''%%0.8g''))');        
            num = num+1;       
    end    
    if isfield(MC.numeric_handles,'SParam2')
            label{num} = sprintf('set(MC.numeric_handles.SParam2,''String'',num2str(SParams(2),''%%0.8g''))');        
            num = num+1;       
    end
    if isfield(MC.numeric_handles,'SParam3')
            label{num} = sprintf('set(MC.numeric_handles.SParam3,''String'',num2str(SParams(3),''%%0.8g''))');        
            num = num+1;       
    end
    if isfield(MC.numeric_handles,'SParam4')
            label{num} = sprintf('set(MC.numeric_handles.SParam4,''String'',num2str(SParams(4),''%%0.8g''))');        
            num = num+1;       
    end    
    if isfield(MC.numeric_handles,'SParam5')
            label{num} = sprintf('set(MC.numeric_handles.SParam5,''String'',num2str(SParams(5),''%%0.8g''))');        
            num = num+1;       
    end
    if isfield(MC.numeric_handles,'SParam6')
            label{num} = sprintf('set(MC.numeric_handles.SParam6,''String'',num2str(SParams(6),''%%0.8g''))');        
            num = num+1;       
    end
    if isfield(MC.numeric_handles,'SParam7')
            label{num} = sprintf('set(MC.numeric_handles.SParam7,''String'',num2str(SParams(7),''%%0.8g''))');        
            num = num+1;       
    end    
    if isfield(MC.numeric_handles,'SParam8')
            label{num} = sprintf('set(MC.numeric_handles.SParam8,''String'',num2str(SParams(8),''%%0.8g''))');        
            num = num+1;       
    end
    if isfield(MC.numeric_handles,'SParam9')
            label{num} = sprintf('set(MC.numeric_handles.SParam9,''String'',num2str(SParams(9),''%%0.8g''))');        
            num = num+1;       
    end
    if isfield(MC.numeric_handles,'SParam10')
            label{num} = sprintf('set(MC.numeric_handles.SParam10,''String'',num2str(SParams(10),''%%0.8g''))');        
            num = num+1;       
    end    
    if isfield(MC.numeric_handles,'SParam11')
            label{num} = sprintf('set(MC.numeric_handles.SParam11,''String'',num2str(SParams(11),''%%0.8g''))');        
            num = num+1;       
    end
    if isfield(MC.numeric_handles,'SParam12')
            label{num} = sprintf('set(MC.numeric_handles.SParam12,''String'',num2str(SParams(12),''%%0.8g''))');        
            num = num+1;       
    end
    

if isfield(MC.numeric_handles,'T')
    label{num}='set(MC.numeric_handles.T,''String'',num2str(T))';    
    num = num +1;    
end

if isfield(MC.numeric_handles,'eps1')
    label{num}='set(MC.numeric_handles.eps1,''String'',num2str(eps1))';    
    num = num +1;    
end

if (gds.options.Userfunctions==1)
    dimu = size(gds.userfunction,2);%       if isfield(MC.numeric_handles,'user_1')    
    if dimu >0 && isfield(MC.numeric_handles,'user')    
        for k = 1:dimu        
            if (gds.options.UserfunctionsInfo(k).state==1)            
                label{num}=sprintf('set(MC.numeric_handles.user_%d,''String'',num2str(hout(%d,i(end)),''%%0.8g''))',k,2+k);                
                num = num+1;               
            end
        end
    end
else dimu = 0;
end

if size(gds.SParams,1)~=0
    if (HTHSNds.index == 1)
        if numsing(1)~=0 && isfield(MC.numeric_handles,'Prod_asp')
            label{num}=sprintf('set(MC.numeric_handles.Prod_asp,''String'',num2str(hout(%d,i(end)),''%%0.8g''))',numsing(1)+dimu);           
            num = num + 1;    
        end
    end
end

if (HTHSNds.index == 2) || (size(gds.SParams,1) == 0)
    if numsing(1)~=0 && isfield(MC.numeric_handles,'disteps1')
        label{num}=sprintf('set(MC.numeric_handles.disteps1,''String'',num2str(hout(%d,i(end)),''%%0.8g''))',numsing(1)+dimu);           
        num = num + 1;    
    end
end


if (gds.options.Eigenvalues==1)&& isfield(MC.numeric_handles,'Re_1')%multipliers numeric window               
    for k = 1:gds.dim    
        label{num}=sprintf('set(MC.numeric_handles.Re_%d,''String'',real(fout(end-gds.dim+%d,i(end))))',k,k);        
        label{num+1}=sprintf('set(MC.numeric_handles.Im_%d,''String'',imag(fout(end-gds.dim+%d,i(end))))',k,k);       
        num = num+2;        
    end
end
if isfield(MC.numeric_handles,'stepsize')
    label{num}='set(MC.numeric_handles.stepsize,''String'',num2str(cds.h,''%0.8g''))';   
    num = num+1;   
end
if isfield(MC.numeric_handles,'npoints')
    label{num}='set(MC.numeric_handles.npoints,''String'',num2str(npoints,''%0.8g''))';    
end

%-------------------------------------------------------------------------    
function draw(varargin)
global gds HTHSNds MC

ndim = size(gds.parameters,1);
if strcmp (varargin{4},'select')
    x = varargin{3};s = varargin{5}; h = varargin{6}; f = varargin{7};    
else
    file = varargin{3}; load(file);  s(1) = []; s(end) = [];    
end
if ~isempty(s)
    sind = cat(1,s.index);    
end
p = size(x,2);
switch varargin{1}
    case 2    
        plot2 = varargin{2};
        if (plot2==[inf inf])                    
            return;
        end
        hold on; d = axis;
        skew = 0.01*[d(2)-d(1) d(4)-d(3)]; 
        watchon;   axis(d);
        plo2{1}  = 'empt'; plo2{2}  = 'empt'; 
        plo2s{1} = 'empt'; plo2s{2} = 'empt'; dimplot = 0;
        for k = 1:2
            if plot2(1,k)<0 
                if plot2(1,k)>=-ndim
                    pl2 = -plot2(1,k);
                    plo2{k} = gds.parameters{pl2,2}*ones(HTHSNds.tps,p); 
                elseif (plot2(1,k)<-ndim) && (plot2(1,k)>=-size(gds.UParams,1)-ndim)
                    pl2 = -plot2(1,k)-ndim;
                    plo2{k} = gds.UParams{pl2,2}*ones(HTHSNds.tps,p);                  
                elseif (plot2(1,k)<-size(gds.UParams,1)-ndim) && (plot2(1,k)>=-size(gds.UParams,1)-size(gds.SParams,1)-ndim)
                    pl2 = -plot2(1,k)-size(gds.UParams,1)-ndim;
                    plo2{k} = gds.SParams{pl2,2}*ones(HTHSNds.tps,p); 
                end   
            end
            if (plot2(1,k)==0)
                plo2{k} = gds.time{1,2}*ones(HTHSNds.tps,p);
            end
                  
            if plot2(1,k) < -(ndim+size(gds.UParams,1)+size(gds.SParams,1)) && plot2(1,k)>-10000*ndim
                if mod(plot2(1,k)+size(gds.UParams,1)+size(gds.SParams,1)+ndim,2)==1
                    plo2{k} = abs(f(HTHSNds.ntst+1+(-plot2(1,k)-(size(gds.UParams,1)+size(gds.SParams,1)+ndim)+1)/2,:));
                else %modulus
                    plo2{k} = angle(f(HTHSNds.ntst+1+(-plot2(1,k)-(size(gds.UParams,1)+size(gds.SParams,1)+ndim)+1)/2,:));
                end
            elseif plot2(1,k) <= -100000*ndim % case Userfunction
                plo2{k}  = ones(HTHSNds.tps,1)*h(2-plot2(1,k)-100000*ndim,:);      
            end
                                                  
            if strcmp(plo2{k},'empt')==1
                if plot2(1,k) <= gds.dim 
                    dimplot = dimplot+1;
                    plo2{k} = x((0:HTHSNds.tps-1)*HTHSNds.nphase+plot2(1,k),:);
                else
                    plo2{k} = ones(HTHSNds.tps,1)*x(plot2(1,k),:);
                end
            end            
        end
        if dimplot == 0
            plo2{1} = plo2{1}(1,:);plo2{2} = plo2{2}(1,:);
        end
        if ~isempty(s)                
            plo2s{1} = plo2{1}(:,sind);
            plo2s{2} = plo2{2}(:,sind);
            plo2{1}(:,sind) = [];
            plo2{2}(:,sind) = [];
        end         
        if strcmp(varargin{4},'redraw')
            if ~strcmp(plo2{1},'empt')&&~strcmp(plo2{2},'empt')
                tmp = plo2{2}';
                if size(tmp,1) == gds.dim
                    for ii = 1:size(tmp,1)
                        line(plo2{1},tmp(ii,:),'LineStyle','-','Color','b');
                    end
                else                    
                    line(plo2{1},plo2{2},'LineStyle','-','Color','b');
                end
            end
            if ~strcmp(plo2s{1},'empt')&&~strcmp(plo2s{2},'empt')
                tmp = plo2s{2}';
                if size(tmp,1) == gds.dim
                    for ii = 1:size(tmp,1)
                        line(plo2s{1},tmp(ii,:),'Linestyle','-','Marker','.','MarkerEdgeColor','r');
                    end
                else                    
                    line(plo2s{1},plo2s{2},'Linestyle','none','Marker','.','MarkerEdgeColor','r');
                end
            end        
            if ~isempty(s)
                text(plo2s{1}(1,:)+skew(1), plo2s{2}(1,:)+skew(2), cat(1,char(s.label)) );            
            end
        else
            plot(plo2s{1},plo2s{2},'co','LineWidth',4);                     
        end
        hold off;
        watchoff;       
    case 3 %3D
        plo = varargin{2};
        if (plo==[inf inf inf]),return;end
        hold on; d = axis;
        skew = 0.01*[d(2)-d(1) d(4)-d(3) d(6)-d(5)];
        watchon;      axis(d);
        plo3{1}  = 'empt'; plo3{2}  = 'empt'; plo3{3}  = 'empt';
        plo3s{1} = 'empt'; plo3s{2} = 'empt'; plo3s{3} = 'empt'; dimplot = 0; 
        for k = 1:3
            if plo(1,k)<0
                if plo(1,k)>=-ndim %case non active parameter
                    pl2 = -plo(1,k);
                    plo3{k} = gds.parameters{pl2,2}*ones(HTHSNds.tps,p); 
                elseif (plo(1,k)<-ndim) && (plo(1,k)>=-size(gds.UParams,1)-ndim)
                    pl2 = -plo(1,k)-ndim;
                    plo3{k} = gds.UParams{pl2,2}*ones(HTHSNds.tps,p);            
                elseif (plo(1,k)<-size(gds.UParams,1)-ndim) && (plo(1,k)>=-size(gds.UParams,1)-size(gds.SParams,1)-ndim)
                    pl2 = -plo(1,k)-size(gds.UParams,1)-ndim;
                    plo3{k} = gds.SParams{pl2,2}*ones(HTHSNds.tps,p);                    
                elseif plo(1,k) <= -100000*ndim % case Userfunction
                    plo3{k}  = ones(HTHSNds.tps,1)*h(2-plo(1,k)-100000*ndim,:);        
                elseif plo(1,k) < -(ndim+size(gds.UParams,1)+size(gds.SParams,1)) % Case multipliers
                    if mod(plo(1,k)+size(gds.UParams,1)+size(gds.SParams,1)+ndim,2)==1
                        plo3{k} = ones(HTHSNds.tps,1)*abs(f((-plo(1,k)-(size(gds.UParams,1)+size(gds.SParams,1)+ndim)+1)/2,:));
                    else %Mod
                        plo3{k} = ones(HTHSNds.tps,1)*angle(f((-plo(1,k)-(size(gds.UParams,1)+size(gds.SParams,1)+ndim))/2,:));
                    end
                end
            end
            if (plo(1,k)==0) %case time
                plo3{k} = gds.time{1,2}*ones(HTHSNds.tps,p);
            end
            if strcmp(plo3{k},'empt')==1
                if plo(1,k)<=gds.dim       %case coordinate
                    dimplot = dimplot+1;
                    plo3{k}=x((0:HTHSNds.tps-1)*HTHSNds.nphase+plo(1,k),:);
                else
                    plo3{k}=x(plo(1,k)*ones(HTHSNds.tps,1),:);  
                end
            end            
        end
        if dimplot == 0
            plo3{1} = plo3{1}(1,:);plo3{2} = plo3{2}(1,:); plo3{3} = plo3{3}(1,:);
        end 
        if ~isempty(s)                
            plo3s{1} = plo3{1}(:,sind);
            plo3s{2} = plo3{2}(:,sind);
            plo3s{3} = plo3{3}(:,sind);
            plo3{1}(:,sind) = [];
            plo3{2}(:,sind) = [];
            plo3{3}(:,sind) = [];
        end
        if strcmp(varargin{4},'redraw')              
            if ~strcmp(plo3{1},'empt')&&~strcmp(plo3{2},'empt')&&~strcmp(plo3{3},'empt')
                line(plo3{1},plo3{2},plo3{3},'linestyle','-','color','b');
            end
            if ~strcmp(plo3s{1},'empt')&&~strcmp(plo3s{2},'empt')&&~strcmp(plo3s{3},'empt')
                plot3(plo3s{1},plo3s{2},plo3s{3},'linestyle','none','Marker','.','MarkerEdgecolor','r');
            end 
            if ~isempty(s)
                text(plo3s{1}(1,:)+skew(1), plo3s{2}(1,:)+skew(2),plo3s{3}(1,:)+skew(3),cat(1,char(s.label)) );            
            end
        else
            plot3(plo3s{1},plo3s{2},plo3s{3},'co','LineWidth',4);                            
        end
        hold off;
        watchoff;
    case 4 %numeric
        parameters = zeros(1,ndim);
        for d=1:ndim
            parameters(d) = gds.parameters{d,2};
        end                
        gds.options.ActiveParams(1,1);
        
        dimu = size(gds.UParams,1);
        UParams = zeros(1,dimu);
        for d=1:dimu
            UParams(d) = gds.UParams{d,2};
        end                
        gds.options.ActiveUParams(1,1);
        
        dims = size(gds.SParams,1);
        SParams = zeros(1,dims);
        for d=1:dims
            SParams(d) = gds.SParams{d,2};
        end                
        gds.options.ActiveSParams(1,1);
        
        T = zeros(1,1);
        T = gds.T;
        gds.options.ActiveT(1,1);
        
        eps1 = zeros(1,1);
        eps1 = gds.eps1;
        gds.options.ActiveEps1(1,1);
        
        dat = get(MC.numeric_fig,'Userdata');
        cds.h='?';
        for k=1:size(dat.label,2)
            eval(dat.label{k});  
        end
end

%-------------------------------------------------------------------------
function output(numsing,xout,s, hout,fout,i)
global gds HTHSNds MC cds

npoints=i(end);
ndim = size(gds.parameters,1);    
if i(1)>3, jj = 1;else jj = i(1)-1;end
p = size(i(1)-jj:i(end),2);
i = i(1)-jj:i(end);
if strcmp(s.label,'00')|strcmp(s.label,'99')
    sind=[];    
else
    sind = find(i==s.index);    
end
if (size(MC.D2,2)>0)  %case 2D
    for siz = 1:size(MC.D2,2)    
        figure(MC.D2(siz));        
        dat = get(MC.D2(siz),'UserData');        
        s1 = eval(dat.label{1});        
        s2 = eval(dat.label{2});        
        if ~isempty(sind)               
            s3 = s1(:,sind);            
            s4 = s2(:,sind);                          
            s1(:,sind) = [];            
            s2(:,sind) = [];            
            line(s1,s2,'Color','b','Linestyle','-');             
            line(s3,s4,'Color','r','Linestyle','-');            
            d = axis;            
            skew = 0.01*[d(2)-d(1) d(4)-d(3)];            
            text(s3(1,:)+skew(1), s4(1,:)+skew(2), s.label );            
        else
            if size(s2,1) == gds.dim            
                for iii = 1:size(s2,1)                
                    line(s1,s2(iii,:),'Color','b','Linestyle','-');                                
                end
            else
                line(s1,s2,'Color','b','Linestyle','-');                
            end
        end
        drawnow;        
    end
end

if (size(MC.D3,2)>0)%case 3D    
    for siz=1:size(MC.D3,2)    
        figure(MC.D3(siz));        
        dat = get(MC.D3(siz),'UserData');          
        s1 = eval(dat.label{1});        
        s2 = eval(dat.label{2});        
        s3 = eval(dat.label{3});        
        if ~isempty(sind)                              
            s4 = s1(:,sind);            
            s5 = s2(:,sind);                                    
            s6 = s3(:,sind);                                   
            s1(:,sind) = [];                                    
            s2(:,sind) = [];            
            s3(:,sind) = [];            
            line(s1,s2,s3,'Color','b','linestyle','-');            
            line(s4,s5,s6,'linestyle','-','color','r');            
            d = axis;                    
            skew = 0.01*[d(2)-d(1) d(4)-d(3) d(6)-d(5)];                    
            text(s4(1,:)+skew(1),s5(1,:)+skew(2), s6(1,:)+skew(3),s.label );            
        else
            line(s1,s2,s3,'linestyle','-','color','b');            
        end
        drawnow;        
    end
end

if ~isempty(MC.numeric_fig) %numeric window is open
    if HTHSNds.index == 1
        parameters = zeros(1,ndim);    
        for d = 1:ndim    
            parameters(d) = gds.parameters{d,2};        
        end

        dimu = size(gds.UParams,1);    
        UParams = zeros(1,dimu);    
        for d = 1:dimu    
            UParams(d) = gds.UParams{d,2};        
        end
        aup = gds.options.ActiveUParams(1,:);    
        for k = 1:length(aup)    
            j = aup(k);                 
            UParams(j) = xout(HTHSNds.ncoords+k,i(end));                        
        end

        dims = size(gds.SParams,1);    
        SParams = zeros(1,dims);    
        for d = 1:dims    
            SParams(d) = gds.SParams{d,2};        
        end
        asp = gds.options.ActiveSParams(1,:);    
        for k = 1:length(asp)    
            j = asp(k);                 
            SParams(j) = xout(HTHSNds.ncoords+length(aup)+k,i(end));                        
        end

        eps1 = xout(HTHSNds.ncoords+length(aup)+length(asp)+1,i(end));
        T = gds.T;   
        eps0 = gds.eps0;
            
    else 
        parameters = zeros(1,ndim);    
        for d = 1:ndim    
            parameters(d) = gds.parameters{d,2};        
        end       
        
        dims = size(gds.SParams,1); 
        SParams = zeros(1,dims);    
        for d = 1:dims
            SParams(d) = gds.SParams{d,2};        
        end        
        
        newind = HTHSNds.ncoords+1;
               
        T = xout(newind,i(end));    
        newind = newind+1;    
        eps0 = gds.eps0;    
        eps1 = xout(newind,i(end));    
        newind = newind+1;
    end
        
    
    dat = get(MC.numeric_fig,'Userdata');       
    for k=1:size(dat.label,2)    
        eval(dat.label{k});          
    end    
end
%-------------------------------------------------------------------------

function load_point(index,x,string,file,varargin)
global gds HTHSNds

load(file);  

if HTHSNds.index == 1
    for i = 1:gds.dim
        gds.coordinates{i,2} = x(i,index);    
    end
    gds.discretization.ntst = HTHSNds.ntst;
    gds.discretization.ncol = HTHSNds.ncol;
    ndim = size(gds.parameters,1);
    for i = 1:ndim
        gds.parameters{i,2} = HTHSNds.P0(i,1);    
    end

    gds.T = HTHSNds.T;
    gds.eps0 = HTHSNds.eps0;
    gds.x0 = HTHSNds.x0;
    
    gds.UParams = []; 
    for i = 1:size(HTHSNds.UParams,1)
        gds.UParams{i,1} = strcat('UParam',num2str(i));        
        gds.UParams{i,2} = HTHSNds.UParams(i);
    end
    for i = 1:size(HTHSNds.ActiveUParams,2)
        gds.UParams{HTHSNds.ActiveUParams(i),2} = x(HTHSNds.ncoords+i,index);        
    end
    gds.SParams = []; 
    for i = 1:size(HTHSNds.SParams,1)
        gds.SParams{i,1} = strcat('SParam',num2str(i));        
        gds.SParams{i,2} = HTHSNds.SParams(i);
    end
    for i = 1:size(HTHSNds.ActiveSParams,2)
        gds.SParams{HTHSNds.ActiveSParams(i),2} = x(HTHSNds.ncoords+size(HTHSNds.ActiveUParams,2)+i,index);               
    end

    gds.eps1 = x(end,index);    
    
            
else
        
    for i = 1:gds.dim
        gds.coordinates{i,2} = x(i,index);   
    end
    gds.discretization.ntst = HTHSNds.ntst;
    gds.discretization.ncol = HTHSNds.ncol;
    
    ndim = size(HTHSNds.P0,1);
    for i = 1:ndim
        gds.parameters{i,2} = HTHSNds.P0(i,1);   
    end
    
    newind = HTHSNds.ncoords+1;
    
    gds.T = x(newind,index);
    newind = newind+1;
    gds.eps1 = x(newind,index);
    gds.period = gds.T*2;
    gds.x0 = HTHSNds.x0;
    gds.eps0 = HTHSNds.eps0;
    
    p= vertcat(gds.parameters{:,2});
    A = cjac(HTHSNds.func,HTHSNds.Jacobian,gds.x0,num2cell(p),[]);
    D = eig(A);

    [V,i] = min(abs(D));
    B = D([1:i-1 i+1:end]);
    nneg = sum(real(B) < 0);
    if (nneg == HTHSNds.nphase-1)
        if min(abs(real(B))) < 1e-10
            nneg = nneg -1;    
        end
    end
    if (nneg == 0)
        if min(abs(real(B))) < 1e-10
            nneg = nneg +1;
        end
    end

    npos = HTHSNds.nphase-nneg-1;   
    
    if nneg == 0
        QS = eye(HTHSNds.nphase);
    else
        [evc0,evl0] = eig(A);
        [K,i] = min(abs(diag(evl0)));
        evc = evc0(:,[1:i-1 i+1:end]);%eigenvectoren zonder eigenvector horende bij 0
        evl = evl0([1:i-1 i+1:end],[1:i-1 i+1:end]);%eigenwaarden zonder 0
    
        val = 0;
        pos = 0;
        k = 1;
        for i=1:size(A,1)-1
            if real(evl(i,i)) < 0
                val(k) = real(evl(i,i));
                pos(k) = i; %geeft je de posities waar de gezochte eigenwaarden en eigenvectoren staan
                k = k+1;
            end
        end
        evcs = evc(:,pos); %de eigenvectoren die meetellen voor de negatieve
        evls = 0;
        for t = 1:size(val,2)
            evls(t) = evl(pos(t),pos(t));
        end

        [a,b] = sort(val); % dan is a(i) = val(b(i))
        VU = zeros(size(A,1),size(val,2));
        sortedevls = zeros(size(pos,2),1);
        for l = 1:size(val,2)
            VU(:,l) = evcs(:,b(l));
            sortedevls(l) = evls(b(l)); %de eigenwaarden gesorteerd van klein naar groot (enkel negatief reeel deel)
        end   

        VU_1 = zeros(size(A,1),size(val,2));
        sortedevls_1 = zeros(size(pos,2),1);
        for f = 1:size(val,2)
            VU_1(:,f) = VU(:,end-f+1);
            sortedevls_1(f) = sortedevls(end-f+1);
        end
        sortedevls = sortedevls_1;

        B = VU_1;
        k = 1;
        while k <= size(B,2)
            ok = 1;
            init = 1;
            while ok == 1 && init <= size(B,1)
                if imag(B(init,k))                     
                    tmp = B(:,k);
                    B(:,k) = real(tmp);
                    B(:,k+1) = imag(tmp);
                    ok = 0;
                end
                init = init+1;
            end
            if ok == 1
               k = k+1; 
            else
                k = k+2;
            end            
        end

        % Compute orthonormal basis for the eigenspace
        [QS,RS] = qr(B);
    end
    
    [V,D] = eig(A);
    [Y,i] = min(abs(diag(D)));
    gds.SParams = [];
    if npos~=0
        Qtot = [QS(:,1:nneg) V(:,i)];
        [Qn,Rn] = qr(Qtot);
        sparams = zeros(1,npos);
        for i=1:npos        
            sparams(i) = 1/gds.eps1*(x(HTHSNds.ncoords-HTHSNds.nphase+1:HTHSNds.ncoords,index)-gds.x0)'*QS(:,nneg+1+i);
        end          

        for i = 1:length(sparams)
            gds.SParams{i,1} = strcat('SParam',num2str(i));
            gds.SParams{i,2} = sparams(i);
        end
    end

    if npos == 0    
        QU = eye(HTHSNds.nphase);                                                      
    else
        [evc0,evl0] = eig(A);
        [K,i] = min(abs(diag(evl0)));
        evc = evc0(:,[1:i-1 i+1:end]);%eigenvectoren zonder eigenvector horende bij 0
        evl = evl0([1:i-1 i+1:end],[1:i-1 i+1:end]);%eigenwaarden zonder 0
                       
        val = 0;                            
        pos = 0;                            
        k = 1;                            
        for i=1:size(A,1)-1                            
            if real(evl(i,i)) > 0                                        
                val(k) = real(evl(i,i));                                                    
                pos(k) = i; %geeft je de posities waar de gezochte eigenwaarden en eigenvectoren staan                                                    
                k = k+1;                                                    
            end
        end
        evcs = evc(:,pos); %de eigenvectoren die meetellen        
        evls = 0;                           
        for t = 1:size(val,2)                            
            evls(t) = evl(pos(t),pos(t));                                        
        end

        [a,b] = sort(val); % dan is a(i) = val(b(i))        
        VU = zeros(size(A,1),size(val,2));                            
        sortedevls = zeros(size(pos,2),1);                           
        for l = 1:size(val,2)                            
            VU(:,l) = evcs(:,b(l));                                        
            sortedevls(l) = evls(b(l)); %de eigenwaarden gesorteerd van klein naar groot (enkel positief reeel deel)                                        
        end

        B = VU;        
        k = 1;                            
        while k <= size(B,2)                            
            ok = 1;                                        
            init = 1;                                       
            while ok == 1 && init <= size(B,1)                                       
                if imag(B(init,k))                                                                         
                    tmp = B(:,k);                                                                
                    B(:,k) = real(tmp);                                                                
                    B(:,k+1) = imag(tmp);                                                                
                    ok = 0;                                                                
                end
                init = init+1;               
            end
            if ok == 1            
                k = k+1;                                                     
            else
                k = k+2;                
            end
        end
        % Compute orthonormal basis for the eigenspace        
        [QU,RU] = qr(B);                            
    end

    [V,D] = eig(A);
    [Y,i] = min(abs(diag(D)));
    gds.UParams = [];
    Qtot = [QU(:,1:npos) V(:,i)];    
    [Qn,Rn] = qr(Qtot);    
    uparams = zeros(1,npos+1);    
    for i = 1:(npos+1)    
        uparams(i) = 1/gds.eps0*(x(1:HTHSNds.nphase,index)-gds.x0)'*Qn(:,i);        
    end
    
    som = sum(uparams.^2);            
    if abs(som-1)>1e-3    
        error('It is not possible to start from this point.')                           
    end

    for i=1:length(uparams)    
        gds.UParams{i,1} = strcat('UParam',num2str(i));        
        gds.UParams{i,2} = uparams(i);        
    end

    if ~isfield(HTHSNds,'TestTolerance') || isempty(HTHSNds.TestTolerance)
        HTHSNds.TestTolerance = 1e-5;
    end

    if npos>1
        vec = [];
        vec = abs(vertcat(gds.SParams{:,2})) < HTHSNds.TestTolerance;    
        aantal = sum(vec);
        if aantal > 0
            HTHSNds.index = 1;        
        end
    end             
    
end

if strcmp(string,'save')
    num = varargin{1};        
    if ~exist('v')
        v = [];        
    end
    if ~exist('s')
        s = [];
    end
    if ~exist('h')
        h = [];
    end
    if ~exist('f')    
        f = [];        
    end
    if ~exist('cds')    
        cds = [];        
    end
    if ~exist('ctype')    
        ctype = [];        
    end
    if ~exist('point')    
        point = [];        
    end
    save(file,'x','v','s','h','f','num','cds','HTHSNds','ctype','point');    
end

%-------------------------------------------------------------------------
function singularities
global gds
  
if size(gds.options.IgnoreSingularity,2)==4
    gds.options = contset(gds.options,'Singularities',0);   
else
    gds.options = contset(gds.options,'Singularities',1);    
end      
%-------------------------------------------------------------------------
    
function start_cont(varargin)
global gds path_sys cds HTHSNds MC

if ~isempty(varargin)
    file = fullfile(path_sys,gds.system,gds.diagram,strcat(gds.curve.new,'.mat'));    
    if exist(file,'file')    
        load(file);        
    else
        errordlg('It is not possible to extend the current curve!','Error extend HTHSN');        
        return         
    end
    [x,v,s,h,f] = cont(x,v,s,h,f,cds);    
else
    if ~isempty(gds.curve.old)    
        file = strcat(gds.curve.old,'.mat');                     
        file = fullfile(path_sys,gds.system,gds.diagram,file);                     
        if exist(file,'file')               
            load(file);             
        end
    end

    systemhandle = str2func(gds.system);    
    
    number = 0; %number of sparams equal to zero
    for i = 1:size(gds.SParams,1)
        if gds.SParams{i,2} == 0
            number = number+1;
        end
    end
   
    if (size(gds.SParams,1)-number) > 0
        n_paru = size(gds.options.ActiveUParams,2);
        n_pars = size(gds.options.ActiveSParams,2);
        dim_npos = size(gds.UParams,1)-1;
        
        for i = 1:n_pars
            if gds.SParams{gds.options.ActiveSParams(i),2} == 0                
                error('A zero SParam should not be denoted as free');
            end        
        end       
        
        if n_paru+n_pars ~= dim_npos+2            
            error('The wrong number of free connection parameters is denoted');
        end   
        
        n_par = size(gds.options.ActiveParams,2);
        if n_par ~= 0
            error('The wrong number of free parameters is denoted');
        end
                                   
        if ~(gds.extravec(3)==1)
            error('eps1 must be free');
        end
        
        if gds.extravec(1)==1
            error('T must not be free');
        end
      
        if HTHSNds.index == 0 
            HTHSNds.index = 1;
            parU = vertcat(gds.UParams{:,2});            
            %parU: alle UParams onder elkaar            
            parS = vertcat(gds.SParams{:,2});                            
            %parS: ALLE SParams onder elkaar                 
            x = [x;gds.x0];     
            par = vertcat(gds.parameters{:,2});
            [x0,v0] = init_HTHSN_HTHSN(systemhandle,x,v,s,par,parU,gds.options.ActiveUParams,parS,gds.options.ActiveSParams,gds.discretization.ntst,gds.discretization.ncol,gds.T,gds.eps1,gds.eps1tol);                                    
        elseif HTHSNds.index == 1
            x0 = x(:,s(num).index);%enkel die éne kolom is belangrijk!                                  
            x0 = [x0(1:(HTHSNds.ntst*HTHSNds.ncol+1)*gds.dim);gds.x0];
            v0 = v(:,s(num).index);            
            parU = vertcat(gds.UParams{:,2});            
            parS = vertcat(gds.SParams{:,2});  
            par = vertcat(gds.parameters{:,2});
            [x0,v0] = init_HTHSN_HTHSN(systemhandle, x0, v0, s(num),par,parU, gds.options.ActiveUParams, parS, gds.options.ActiveSParams, gds.discretization.ntst,gds.discretization.ncol,gds.T,gds.eps1,gds.eps1tol);            
        end   
        
    else
        
        if ~(size(gds.options.ActiveUParams,2) == 0)
            error('The wrong number of free unstable parameters is denoted');
        end
        
        if ~(size(gds.options.ActiveSParams,2) == 0)
            error('The wrong number of free stable parameters is denoted');
        end
        
        n_par = size(gds.options.ActiveParams,2);
        if n_par ~= 0
            error('0 free system parameters are needed');
        end      
        
        if ~(gds.extravec(1)==1)
            error('T must be free');
        end       
        if ~(gds.extravec(3)==1)
            error('eps1 must be free');
        end        
        
        if HTHSNds.index == 0
            HTHSNds.index = 2;     
            x = [x;gds.x0];                        
            par = vertcat(gds.parameters{:,2});
            [x0,v0] = init_HTHSN_HTHSN(systemhandle,x,v,s,par,[],[],[],[],gds.discretization.ntst,gds.discretization.ncol,gds.T,gds.eps1,gds.eps1tol);         
        elseif HTHSNds.index == 1
            HTHSNds.index = 2;
            x0 = x(:,s(num).index);            
            v0 = v(:,s(num).index);                
            x0 = [x0(1:(HTHSNds.ntst*HTHSNds.ncol+1)*gds.dim);gds.x0];
            par = vertcat(gds.parameters{:,2});            
            gds.T = gds.period/2;            
            [x0,v0] = init_HTHSN_HTHSN(systemhandle, x0, v0, s(num), par,[],[],[],[],gds.discretization.ntst,gds.discretization.ncol,gds.T,gds.eps1,gds.eps1tol);
        else
            x0 = x(:,s(num).index);            
            v0 = v(:,s(num).index);            
            x0 = [x0(1:(HTHSNds.ntst*HTHSNds.ncol+1)*gds.dim);gds.x0];
            par = vertcat(gds.parameters{:,2});            
            gds.T = gds.period/2;            
            [x0,v0] = init_HTHSN_HTHSN(systemhandle, x0, v0, s(num), par,[],[],[],[],gds.discretization.ntst,gds.discretization.ncol,gds.T,gds.eps1,gds.eps1tol);
        end
        
    end                       
    if isempty(x0)    
        return        
    end
    [x,v,s,h,f] = cont(@homotopysaddlenode,x0,v0,gds.options);    
end
file = fullfile(path_sys,gds.system,gds.diagram,gds.curve.new);
point = gds.point; ctype = gds.type;
num = 1;
%     gds.curve.old = gds.curve.new;  door dit weg te laten kan je weer
%     beginnen vanaf die tempadwgyk_cyclecon, anders moet je telkens
%     opnieuw de tijdsintegratie doen!!!
status = mkdir(path_sys,gds.system);
dir = fullfile(path_sys,gds.system);
status = mkdir(dir,gds.diagram);
if ~isempty(x)
    save(file,'x','v','s','h','f','cds','HTHSNds','ctype','point','num');    
end    

%-------------------------------------------------------------------------

function j = start_USparams(handles,j)
global gds HTHSNds
color = [1 1 1];

j = j-2;
pos = [5 j 38 1.80];user.num = 0; user.pos = pos;
string = 'Connection parameters';
stat = uicontrol(handles.figuur,'Style','text','String',string,'Tag','Connection parameters','BackGroundColor',color,'Units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',pos,'UserData',user);

dim_up = size(gds.UParams,1);
for i = 1:dim_up
    %eerst voor de c's
    if strcmp(gds.UParams{i,1},'')
        string = '';
    else
        string = num2str(gds.UParams{i,2},'%0.8g');
    end
    value = 0;
    if ~isempty(gds.options.ActiveUParams)
        str = gds.options.ActiveUParams;
        x = str(str==i);
        if ~isempty(x)
            value = 1;
        end
    end
    
    tag2 = strcat('edit',num2str(gds.dim+size(gds.parameters,1)+i));
    user.num = num2str(i);
    edit = uicontrol(handles.figuur,'Style','edit','HorizontalAlignment','left','String',string,'Tag',tag2,'Backgroundcolor',color,'units','characters','fontname','FixedWidth','fontsize',12);
    rad = uicontrol(handles.figuur,'Style','radiobutton','String',gds.UParams{i,1},'Tag','ActiveUParams','Callback','set_option','Max',1,'Value',value,'BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
    set(edit,'Callback','HT_Callback');
    j = j-2;
    pos = [20 j 25 1.8]; user.pos=pos;
    set(edit,'Position',pos,'UserData',user);
        
    pos = [2 j 18 1.8];user.pos=pos;
    set(rad,'Position',pos,'UserData',user);
end

dim_nneg = size(gds.SParams,1);

if ~isfield(HTHSNds,'TestTolerance') || isempty(HTHSNds.TestTolerance)
    HTHSNds.TestTolerance = 1e-5;
end

for i = 1:dim_nneg
    prod = 1;
    stablepars = [];
    index = [];
    for i = 1:dim_nneg
        if abs(gds.SParams{i,2})~=0
            prod = prod*abs(gds.SParams{i,2});
            stablepars = [stablepars;abs(gds.SParams{i,2})];
            index = [index;i];
        end
    end
    if prod < HTHSNds.TestTolerance
        [val,pos] = min(stablepars);
        gds.SParams{index(pos),2} = 0;
    end
end


for i = 1:dim_nneg
    if abs(gds.SParams{i,2}) < HTHSNds.TestTolerance
        gds.SParams{i,2} = 0;
    end
end

for i = 1:dim_nneg
    %nu voor de tau's
    
    if strcmp(gds.SParams{i,1},'')
        string = '';
    else
        string = num2str(gds.SParams{i,2},'%0.8g');
    end
    value = 0;
    if ~isempty(gds.options.ActiveSParams)
        str = gds.options.ActiveSParams;
        x = str(str==i);
        if ~isempty(x)
            value = 1;
        end
    end
    
    tag2 = strcat('edit',num2str(gds.dim+size(gds.parameters,1)+dim_up+i));
    user.num = num2str(i);
    edit = uicontrol(handles.figuur,'Style','edit','HorizontalAlignment','left','String',string,'Tag',tag2,'Backgroundcolor',color,'units','characters','fontname','FixedWidth','fontsize',12);
    rad = uicontrol(handles.figuur,'Style','radiobutton','String',gds.SParams{i,1},'Tag','ActiveSParams','Callback','set_option','Max',1,'Value',value,'BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
    set(edit,'Callback','HT_Callback');
    j = j-2;
    pos = [20 j 25 1.8]; user.pos=pos;
    set(edit,'Position',pos,'UserData',user);
        
    pos = [2 j 18 1.8];user.pos=pos;
    set(rad,'Position',pos,'UserData',user);
    
end
    
guidata(handles.figuur,handles);

%--------------------------------------------------------------------------
function j = start_eps1tol(handles,j)
global gds
color = [1 1 1];

j = j-2;
pos = [5 j 38 1.80];user.num = 0; user.pos = pos;
string = 'eps1 tolerance';
stat = uicontrol(handles.figuur,'Style','text','String',string,'Tag','eps1 tolerance','BackGroundColor',color,'Units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',pos,'UserData',user);

tag  = strcat('edit',num2str(1001));
%string = num2str(gds.T);
if ~isfield(gds,'eps1tol') || isempty(gds.eps1tol)
    gds.eps1tol = 1e-2;
end
value = gds.eps1tol;
edit  = uicontrol(handles.figuur,'Style','edit','HorizontalAlignment','left','String',gds.eps1tol,'Tag',tag,'BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
stat = uicontrol(handles.figuur,'Style','text','HorizontalAlignment','left','String','eps1tol','Tag','eps1tol','Callback','set_option','Max',1,'Value',value,'Backgroundcolor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(edit,'Callback','HT_Callback');
j = j-2;
pos = [20 j 25 1.8]; user.num = 0; user.pos = pos;
set(edit,'Position',pos,'UserData',user);

pos = [2 j 18 1.8];user.pos = pos; 
set(stat,'Position',pos,'UserData',user);
guidata(handles.figuur,handles);


%--------------------------------------------------------------------------
function j = start_period(handles,j)
global gds;
color = [1 1 1];

j = j-2;
pos = [5 j 38 1.80];user.num = 0; user.pos = pos;
string = 'Homoclinic parameters';
stat = uicontrol(handles.figuur,'Style','text','String',string,'Tag','Homoclinic parameters','BackGroundColor',color,'Units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',pos,'UserData',user);


tag  = strcat('edit',num2str(999));
%string = num2str(gds.T);
value = 0;
edit  = uicontrol(handles.figuur,'Style','edit','HorizontalAlignment','left','String',gds.T,'Tag',tag,'BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
rad = uicontrol(handles.figuur,'Style','radiobutton','HorizontalAlignment','left','String','T','Tag','Tpar','Callback','set_option','Max',1,'Value',value,'Backgroundcolor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(edit,'Callback','HT_Callback');
j = j-2;
pos = [20 j 25 1.8]; user.num = 0; user.pos = pos;
set(edit,'Position',pos,'UserData',user);

pos = [2 j 18 1.8];user.pos = pos; 
set(rad,'Position',pos,'UserData',user);
guidata(handles.figuur,handles);

tag  = strcat('edit',num2str(1000));
string = num2str(gds.eps1);
edit = uicontrol(handles.figuur,'Style','edit','HorizontalAlignment','left','String',string,'Tag',tag,'Backgroundcolor',color,'units','characters','fontname','FixedWidth','fontsize',12);
rad  = uicontrol(handles.figuur,'Style','radiobutton','HorizontalAlignment','left','String','eps1','Tag','eps1par','Callback','set_option','Max',1,'BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(edit,'Callback','HT_Callback');

j = j-2;
pos = [20 j 25 1.8];user.pos = pos; 
set(edit,'Position',pos,'UserData',user);
pos = [2 j 18 1.8]; user.pos = pos;
set(rad,'Position',pos,'UserData',user);
guidata(handles.figuur,handles);

%--------------------------------------------------------------------------
function j=start_testfunctions(handles,j)
global gds HTHSNds

if ~isfield(HTHSNds,'index') || isempty(HTHSNds.index)
    
else
    if ~isempty(gds.SParams)      
        if (HTHSNds.index == 0 || HTHSNds.index == 1)        
            sp = vertcat(gds.SParams{:,2});
            if ~isfield(HTHSNds,'TestTolerance') || isempty(HTHSNds.TestTolerance)                
                vec = zeros(size(sp))+1e-8;
            else
                vec = zeros(size(sp))+HTHSNds.TestTolerance;
            end

            if abs(sp) < vec
            else
                color = [1 1 1];
                j = j-2;
                pos = [5 j 38 1.80];user.num = 0; user.pos = pos;
                string = strcat('Testfunctions (',gds.type,')');
                stat = uicontrol(handles.numeric_fig,'Style','text','String',string,'Tag','testfunctions','BackGroundColor',color,'Units','characters','fontname','FixedWidth','fontsize',12);
                set(stat,'Position',pos,'UserData',user);

                j = j-2;
                post = [2 j 20 1.8]; user.pos = post;
                pose = [20 j 24 1.8];            
                tag = 'Prod_asp';
                stat = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String',tag,'BackGroundColor',color,'Units','characters','fontname','FixedWidth','fontsize',12);
                edit = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','','Tag',tag,'BackGroundColor',color,'Units','characters','fontname','FixedWidth','fontsize',12);
                set(stat,'Position',post,'UserData',user); user.pos = pose;
                set(edit,'Position',pose,'UserData',user);
            end
        end
    end
    if (HTHSNds.index == 2) || (size(gds.SParams,1)==0)
        color = [1 1 1];      
        j = j-2;    

        pos = [5 j 38 1.80];user.num = 0; user.pos = pos;    
        string = strcat('Testfunctions (',gds.type,')');    
        stat = uicontrol(handles.numeric_fig,'Style','text','String',string,'Tag','testfunctions','BackGroundColor',color,'Units','characters','fontname','FixedWidth','fontsize',12);    
        set(stat,'Position',pos,'UserData',user);    
        j = j-2;    
        post = [2 j 20 1.8]; user.pos = post;    
        pose = [20 j 24 1.8];                
        tag = 'disteps1';    
        stat = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','eps1-eps1tol','BackGroundColor',color,'Units','characters','fontname','FixedWidth','fontsize',12);    
        edit = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','','Tag',tag,'BackGroundColor',color,'Units','characters','fontname','FixedWidth','fontsize',12);    
        set(stat,'Position',post,'UserData',user); user.pos = pose;    
        set(edit,'Position',pose,'UserData',user);        
    end
end    

guidata(handles.numeric_fig,handles);

%--------------------------------------------------------------------------
function j = numeric_HTpars(handles,j)
global gds
color = [1 1 1];
j = j-2;
stat = uicontrol(handles.numeric_fig,'Style','text','String','Connection parameters','BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
pos = [5 j 38 1.8];user.num=0;user.pos=pos;
set(stat,'Position',pos,'UserData',user);

dim_s = size(gds.SParams,1);
for i = 1:dim_s
    string = '';
    tag1 = strcat('text_SParam',num2str(i));
    tag2 = strcat('SParam',num2str(i));    
    stat1 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String',gds.SParams{i,1},'Tag',tag1,'BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
    j = j-2;
    pos = [2 j 18 1.8];user.pos=pos;    
    set(stat1,'Position',pos,'UserData',user);

    stat2 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String',string,'Tag',tag2,'BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
    pos = [20 j 24 1.8];user.pos = pos;
    set(stat2,'Position',pos,'UserData',user);
end
guidata(handles.numeric_fig,handles)

%--------------------------------------------------------------------------
function j = numeric_HTeps1(handles,j)
global gds
color = [1 1 1];
j = j-2;
stat = uicontrol(handles.numeric_fig,'Style','text','String','Homoclinic parameters','BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
pos = [5 j 38 1.8];user.num=0;user.pos=pos;
set(stat,'Position',pos,'UserData',user);
j = j-2;
post = [2 j 18 1.8];user.pos=post;
pose = [20 j 21 1.8];
stat1 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','T','Tag','text_T','BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
set(stat1,'Position',post,'UserData',user);user.pos=pose;
stat2 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','','Tag','T','Backgroundcolor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
set(stat2,'Position',pose,'UserData',user);
guidata(handles.numeric_fig,handles);
j = j-2;
post = [2 j 18 1.8];user.pos=post;
pose = [20 j 21 1.8];
stat1 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','eps1','Tag','text_eps1','BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
set(stat1,'Position',post,'UserData',user);user.pos=pose;
stat2 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','','Tag','eps1','Backgroundcolor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
set(stat2,'Position',pose,'UserData',user);
guidata(handles.numeric_fig,handles);
