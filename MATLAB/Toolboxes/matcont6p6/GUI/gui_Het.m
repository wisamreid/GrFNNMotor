function gui_Het
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
    [point,str] = strtok(tag,'_');
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
    else curvetype;
    end
%-------------------------------------------------------------------------    
function curvetype
global gds path_sys MC
    if isempty(gds.point)
        gds.point = 'Het ';
    end
    gds.type = 'Het ';
    set(0,'ShowHiddenHandles','on');        
    h = findobj('Type','uimenu','Tag','forward'); set(h,'enable','on');
    h = findobj('Type','uimenu','Tag','backward'); set(h,'enable','on');
    set(0,'ShowHiddenHandles','off');  
    matcont('make_curve');
    load_matcont;
    file = strcat(gds.system,'.mat');
    file = fullfile(path_sys,file);
    save(file,'gds');
    starter; continuer;
    if  ~isempty(MC.numeric_fig),numeric;end
    set(0,'ShowHiddenHandles','on');        
    h = findobj('Type','uimenu','Tag','window');set(h,'enable','on');
    h = findobj('Type','uimenu','Tag','extend'); set(h,'enable','off');
    set(0,'ShowHiddenHandles','off');        
%-------------------------------------------------------------------------    
function starter1(handles)
global gds hetds
    ndim = size(gds.parameters,1);
    color = [1 1 1];
    gds.options.ActiveUParams = [];
    gds.options.ActiveSParams = [];
    
    s = ndim*2+20+2*7+2;   
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
    pos1 = start_period(handles,pos1);
    pos1 = starter('start_jacobian',handles,pos1);
    pos1 = starter('start_discret',handles,pos1);
    gds.options.IgnoreSingularity = 1:15;
%    pos1 = start_monitor(handles,pos1);
    if isfield(gds,'userfunction')
        starter('start_userfunctions',handles,pos1);
    end
    starter('in_start',handles);
%-------------------------------------------------------------------------    
function load_layout(handles)
global gds 
    if size(gds.numeric.Het,1)<7
         gds.numeric.Het{6,1} = 'user functions';
         gds.numeric.Het{6,2} = 0;
         gds.numeric.Het{7,1} = 'npoints';
         gds.numeric.Het{7,2} = 0;
     end
     for i = 1:7
         if gds.numeric.Het{i,2}==1
             gds.numeric.Het{i,1} = upper(gds.numeric.Het{i,1});
         elseif gds.numeric.Het{i,2}==0
             gds.numeric.Het{i,1} = lower(gds.numeric.Het{i,1});
         end
        string{i,1} = gds.numeric.Het{i,1};
    end    
    set(handles.layoutbox,'String',string);
%-------------------------------------------------------------------------    
function layoutbox(list,index_selected)
global gds 
    c = gds.numeric.Het{index_selected,2};    
    gds.numeric.Het{index_selected,2} = 1-c;    
%-------------------------------------------------------------------------    
function numeric1(handles)
global gds
    ndim = size(gds.parameters,1);
    if ~isfield(gds.numeric,'Het') || strcmp(gds.numeric.Het{4,1},'multipliers') || strcmp(gds.numeric.Het{4,1},'MULTIPLIERS')   
        gds.numeric.Het = {'parameters' 1;'period' 1;'testfunctions' 0;'eigenvalues' 1;'current stepsize' 0};
    end
     if size(gds.numeric.Het,1)<7
         gds.numeric.Het{6,1} = 'user functions';
         gds.numeric.Het{6,2} = 0;
         gds.numeric.Het{7,1} = 'npoints';
         gds.numeric.Het{7,2} = 0;
     end
         gds.numeric.Het{8,1} = 'eps0';
         gds.numeric.Het{8,2} = 1;
         gds.numeric.Het{9,1} = 'eps1';
         gds.numeric.Het{9,2} = 1;
    if isfield(gds,'userfunction')
        dimu = size(gds.userfunction,2);
    else dimu=0; 
    end
    s = 2;
    if gds.numeric.Het{1,2}==1
        s = s+2*ndim+2;
    end
    if gds.numeric.Het{2,2}==1
        s = s+4;
    end
    if gds.numeric.Het{3,2}==1
        s = s+8+2;
    end  
    if gds.numeric.Het{4,2}==1
        s = s+8*gds.dim+2;
    end  
    if gds.numeric.Het{5,2}==1
        s = s+2;
    end 
    if gds.numeric.Het{6,2}==1
        s = s+2*dimu+2;
    end
    if gds.numeric.Het{7,2}==1
        s = s+2;
    end    
    if gds.numeric.Het{8,2}==1
        s = s+2;
    end    
    if gds.numeric.Het{9,2}==1
        s = s+2;
    end    
    slider_step(1) = 2/s;
    slider_step(2) = 2/s;
    set(handles.slider1,'sliderstep',slider_step,'max',s,'min',0,'Value',s);
    j = s;    
    if gds.numeric.Het{1,2}==1
        j = numeric('start_parameters',handles,j);
    end
    if gds.numeric.Het{2,2}==1
       j = numeric_hetpars(handles,j);
    end
    if gds.numeric.Het{3,2}==1
%        j = start_testfunctions(handles,j);
    end
    if gds.numeric.Het{4,2}==1
       j = start_eigenvalues(handles,j);
    end  
    if gds.numeric.Het{5,2}==1
        j = numeric('start_stepsize',handles,j);
    end 
    if gds.numeric.Het{6,2}==1
        j = numeric('start_userfunctions',handles,j);
    end
    if gds.numeric.Het{7,2}==1
        numeric('start_npoints',handles,j);
    end

    numeric('in_numeric',handles);
%-------------------------------------------------------------------------    
function [string,plotopts] = load_draw
global gds
    plotopts = {};
    string = {'Coordinates';'Parameters';'Period'};
    if gds.options.Eigenvalues==1
        string{4} = 'Eigenvalues';
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
global gds hetds
    ndim = size(gds.parameters,1);
    len  = size(pl,2);
    leng = size(x,1);
    e(1,len) = 0;
    
    for j = 1:len
        switch pl(j).type
        case 'Coordinates'
            e(1,j) = pl(j).val;
        case 'Parameters'
            val = 0;
            for k = 1:length(gds.options.ActiveParams)  
                if pl(j).val==gds.options.ActiveParams(1,k)
                    e(1,j) = hetds.PeriodIdx+k;
                    val = 1;
                end
            end
            if val == 0
                e(1,j) = -pl(j).val;
            end
        case 'Eigenvalues'
            if gds.options.Eigenvalues==0
                errordlg('It isn''t possible to plot the eigenvalues(they are not being computed)!');
            else
                e(1,j) = -ndim-pl(j).val;
            end
        case 'Period'
            e(1,j) = hetds.PeriodIdx+length(gds.options.ActiveParams)+1;
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
    label = cell(len,1);
    lab = cell(len,1);
    for k = 1:len
        label{k}='empty';
        lab{k}='empty';
        if plot2(1,k)<0 && plot2(1,k)>=-ndim
            pl2  = -plot2(1,k);
            label{k}  = sprintf('gds.parameters{%d,2}*ones(hetds.tps,p)',pl2);  
            lab{k}  = sprintf('gds.parameters{%d,2}',pl2);      
        end     
        if plot2(1,k)==0
            label{k}  = 'gds.time{1,2}*ones(hetds.tps,p)';
            lab{k} = 'gds.time{1,2}';
        end
        if plot2(1,k) < -ndim %case eigenvalues
            if plot2(1,k) <= -100000*ndim % case Userfunction
                label{k}  = sprintf('hout(2+%d,i)',-plot2(1,k)-100000*ndim); 
                lab{k}  = sprintf('hout(2+%d,i)',-plot2(1,k)-100000*ndim); 
            elseif mod(plot2(1,k)+ndim,2)==0
                label{k}  = sprintf('angle(fout(hetds.ntst+1+%d,i))',(-plot2(1,k)-ndim)/2);            
                lab{k} = sprintf('angle(fout(hetds.ntst+1+%d,i))',(-plot2(1,k)-ndim)/2);            
            else %Re
                label{k} = sprintf('abs(fout(hetds.ntst+1+%d,i))',(-plot2(1,k)-ndim+1)/2);            
                lab{k} = sprintf('abs(fout(hetds.ntst+1+%d,i))',(-plot2(1,k)-ndim+1)/2);                    
            end
        end
        if strcmp(label{k},'empty')==1
            dimplot = dimplot+1;
            if plot2(1,k) <= gds.dim %case coordinate
                label{k}  = sprintf('xout((0:hetds.tps-1)*hetds.nphase+%d,i)',plot2(1,k));            
                lab{k}  = sprintf('xout((0:hetds.tps-1)*hetds.nphase+%d,i)]',plot2(1,k));                      
            else %case active parameter
                label{k} = sprintf('ones(hetds.tps,1)*xout(%d,i)',plot2(1,k));            
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
global gds  MC
    ndim = size(gds.parameters,1);
    num = 1;
    label{ndim} = '';
    if isfield(MC.numeric_handles,'param1')
        for d = 1:ndim
            label{num} = sprintf('set(MC.numeric_handles.param%d,''String'',num2str(parameters(%d),''%%0.8g''))',d,d);
            num = num +1;
        end            
    end
    if isfield(MC.numeric_handles,'T')
        label{num}='set(MC.numeric_handles.T,''String'',num2str(T))';
        num = num +1;
    end
    if isfield(MC.numeric_handles,'eps0')
        label{num}='set(MC.numeric_handles.eps0,''String'',num2str(eps0))';
        num = num +1;
    end
    if isfield(MC.numeric_handles,'eps1')
        label{num}='set(MC.numeric_handles.eps1,''String'',num2str(eps1))';
        num = num +1;
    end
    if (gds.options.Userfunctions==1)
        dimu = size(gds.userfunction,2);
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
 
    if (gds.options.Eigenvalues==1)&& isfield(MC.numeric_handles,'Re_1')%multipliers numeric window
        for k = 1:2*gds.dim
            label{num}=sprintf('set(MC.numeric_handles.Re_%d,''String'',real(fout(end-2*gds.dim+%d,i(end))))',k,k);
            label{num+1}=sprintf('set(MC.numeric_handles.Im_%d,''String'',imag(fout(end-2*gds.dim+%d,i(end))))',k,k);
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
global gds hetds MC
    ndim = size(gds.parameters,1);
     if strcmp (varargin{4},'select')
        x = varargin{3};s = varargin{5}; f = varargin{7};
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
            if plot2(1,k)<0 && plot2(1,k)>=-ndim
                pl2 = -plot2(1,k);
                plo2{k} = gds.parameters{pl2,2}*ones(hetds.tps,p) ; 
            end   
            if (plot2(1,k)==0)
                plo2{k} = gds.time{1,2}*ones(hetds.tps,p);
            end
            if plot2(1,k) <= -100000*ndim % case Userfunction
                plo2{k}  = ones(hetds.tps,1)*h(2-plot2(1,k)-100000*ndim,:);
            
            elseif plot2(1,k) < -ndim
                if mod(plot2(1,k)+ndim,2)==1
                    plo2{k} = abs(f(hetds.ntst+1+(-plot2(1,k)-ndim+1)/2,:));
                else %modulus
                    plo2{k} = angle(f(hetds.ntst+1+(-plot2(1,k)-ndim+1)/2,:));
                end
            end
            if strcmp(plo2{k},'empt')==1
                if plot2(1,k) <= gds.dim 
                    dimplot = dimplot+1;
                    plo2{k} = x((0:hetds.tps-1)*hetds.nphase+plot2(1,k),:);
                else
                    plo2{k} = ones(hetds.tps,1)*x(plot2(1,k),:);
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
            if plo(1,k)<0 && plo(1,k)>=-ndim %case non active parameter
                pl2 = -plo(1,k);
                plo3{k} = gds.parameters{pl2,2}*ones(hetds.tps,p); 
            elseif plo(1,k) <= -100000*ndim % case Userfunction
                plo3{k}  = ones(hetds.tps,1)*h(2-plo(1,k)-100000*ndim,:);                      
            elseif plo(1,k) < -ndim % Case multipliers
                if mod(plo(1,k)+ndim,2)==1
                    plo3{k} = ones(hetds.tps,1)*abs(f((-plo(1,k)-ndim+1)/2,:));
                else %Mod
                    plo3{k} = ones(hetds.tps,1)*angle(f((-plo(1,k)-ndim)/2,:));
                end
            end     
            if (plo(1,k)==0) %case time
                plo3{k} = gds.time{1,2}*ones(hetds.tps,p);
            end
            if strcmp(plo3{k},'empt')==1
                if plo(1,k)<=gds.dim       %case coordinate
                    dimplot = dimplot+1;
                    plo3{k}=x((0:hetds.tps-1)*hetds.nphase+plo(1,k),:);
                else %case active parameter
                    plo3{k}=x(plo(1,k)*ones(hetds.tps,1),:);  
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
        
        T = zeros(1,1);
        T = gds.T;
        gds.options.ActiveT(1,1);
        
        eps0 = zeros(1,1);
        eps0 = gds.eps0;
        gds.options.ActiveEps0(1,1);
        
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
global gds hetds MC cds
    gds.extravec = hetds.extravec;
    ndim = size(gds.parameters,1);
    npoints=i(end);
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
        parameters = zeros(1,ndim);
        for d = 1:ndim
            parameters(d) = gds.parameters{d,2};
        end   
        if ~isempty(gds.options.ActiveParams)
            ap = gds.options.ActiveParams(1,:);
            for k = 1:length(ap)
                parameters(ap(k)) = xout(hetds.PeriodIdx+k,i(end));                
            end        
        end
        newind = hetds.PeriodIdx+length(gds.options.ActiveParams)+1;
        
        if gds.extravec(1)
            T = xout(newind,i(end));
            newind = newind+1;
        else 
            T = gds.T;
        end
        
        if gds.extravec(2)
            eps0 = xout(newind,i(end));
            newind = newind+1;
        else 
            eps0 = gds.eps0;
        end
        
        if gds.extravec(3)
            eps1 = xout(newind,i(end));
            newind = newind+1;
        else
            eps1 = gds.eps1;
        end
        dat = get(MC.numeric_fig,'Userdata');
        for k=1:size(dat.label,2)
            eval(dat.label{k});  
        end
    end
%-------------------------------------------------------------------------    
function load_point(index,x,string,file,varargin)
global gds hetds 
    load(file);
    for i = 1:gds.dim
        gds.coordinates{i,2} = x(i,index);
    end
    gds.discretization.ntst = hetds.ntst;
    gds.discretization.ncol = hetds.ncol;
    ndim = size(hetds.P0,1);
    for i = 1:ndim
        gds.parameters{i,2} = hetds.P0(i,1);
    end    

    for i = 1:length(hetds.ActiveParams)
        gds.parameters{hetds.ActiveParams(i),2} = x(hetds.PeriodIdx+i,index);
    end    
    newind = hetds.PeriodIdx+length(hetds.ActiveParams)+1;
    gds.T = hetds.T;
    if gds.extravec(1)
        gds.T = x(newind,index);
        newind = newind+1;
    end
    gds.period = 2*gds.T;
    gds.eps0 = hetds.eps0;
    if gds.extravec(2)
        gds.eps0 = x(newind,index);
        newind = newind+1;
    end
    gds.eps1 = hetds.eps1;
    if gds.extravec(3)
        gds.eps1 = x(newind,index);
    end
    
    if strcmp(string,'save')
        num = varargin{1};
        save(file,'x','v','s','h','f','num','cds','hetds','ctype','point');
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
global gds path_sys cds hetds MC
     if ~isempty(varargin)
         file = fullfile(path_sys,gds.system,gds.diagram,strcat(gds.curve.new,'.mat'));
         if exist(file,'file')
             load(file);
         else
             errordlg('It is not possible to extend the current curve!','Error extend Het');
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
         switch gds.point           
             case {'HTHet '}
                 par = vertcat(gds.parameters{:,2});
                 x0 = x(:,s(num).index); 
                 x0 = [x0(1:(gds.discretization.ntst*gds.discretization.ncol+1)*gds.dim);gds.x0;gds.x1];  
                 v0 = v(:,s(num).index);
                 [x0,v0] = init_HTHet_Het(systemhandle, x0, v0, s(num), par, gds.options.ActiveParams, gds.discretization.ntst, gds.discretization.ncol,gds.extravec,gds.T,gds.eps0,gds.eps1);
             otherwise %case heteroclinic starts from another heteroclinic                     
                 par = vertcat(gds.parameters{:,2});   
                 x0 = x(:,s(num).index);            
                 x0 = [x0(1:(gds.discretization.ntst*gds.discretization.ncol+1)*gds.dim);gds.x0;gds.x1];  
                 v0 = v(:,s(num).index);
                 [x0,v0] = init_Het_Het(systemhandle,x0,v0,s(num),par,gds.options.ActiveParams,gds.discretization.ntst,gds.discretization.ncol,gds.extravec,gds.T,gds.eps0,gds.eps1);                                                 
         end         
         if isempty(x0)
             return
         end
         [x,v,s,h,f] = cont(@heteroclinic,x0,v0,gds.options);   
     end
     file = fullfile(path_sys,gds.system,gds.diagram,gds.curve.new);
     point = gds.point; ctype = gds.type;     
     num = 1;
     status = mkdir(path_sys,gds.system);
     dir = fullfile(path_sys,gds.system);
     status = mkdir(dir,gds.diagram);
    if ~isempty(x)
        save(file,'x','v','s','h','f','cds','hetds','ctype','point','num');
    end


%---------------------------------------------------------------------------
function j = start_monitor(handles,j)
global gds
color = [1 1 1];
j = j-2;
pos = [10 j 35 1.80]; user.num = 0; user.pos = pos;
stat = uicontrol(handles.figuur,'Style','text','String','Monitor Singularities','Tag','singularities','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',pos,'UserData',user);

guidata(handles.figuur,handles);


%-----------------------------------------------------------------------------------
function j=start_testfunctions(handles,j)
global gds 
color = [1 1 1];
j = j-2;
pos = [5 j 38 1.80];user.num = 0; user.pos = pos;
string = strcat('Testfunctions (',gds.type,')');
stat = uicontrol(handles.numeric_fig,'Style','text','String',string,'Tag','testfunctions','BackGroundColor',color,'Units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',pos,'UserData',user);

guidata(handles.numeric_fig,handles);


%--------------------------------------------------------------------------
function j = start_period(handles,j)
global gds;
color = [1 1 1];

j = j-2;
pos = [5 j 38 1.80];user.num = 0; user.pos = pos;
string = 'Heteroclinic parameters';
stat = uicontrol(handles.figuur,'Style','text','String',string,'Tag','Heteroclinic parameters','BackGroundColor',color,'Units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',pos,'UserData',user);

if ~isfield(gds,'extravec')
    gds.extravec = [1 0 0];
end

j = j-2;
post = [2 j 18 1.8]; user.num = 0; user.pos = post; pose = [20 j 25 1.8];
gds.T = gds.period/2;
tag  = strcat('edit',num2str(1));
edit = uicontrol(handles.figuur,'Style','edit','HorizontalAlignment','left','String',gds.T,'Tag',tag,'Backgroundcolor',color,'units','characters','fontname','FixedWidth','fontsize',12);
value = gds.extravec(1);
rad  = uicontrol(handles.figuur,'Style','radiobutton','String','T','Tag','Tpar','Callback','set_option','Max',1,'Value',value,'BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(edit,'Callback','Hom_Callback');
user.pos = pose; set(edit,'Position',pose,'UserData',user);
    pos = [2 j 18 1.8]; user.pos=pos;
    set(rad,'Position',pos,'UserData',user); guidata(handles.figuur,handles);

j = j-2;
post = [2 j 18 1.8]; user.num = 0; user.pos = post; pose = [20 j 25 1.8];
tag  = strcat('edit',num2str(2));
edit = uicontrol(handles.figuur,'Style','edit','HorizontalAlignment','left','String',gds.eps0,'Tag',tag,'Backgroundcolor',color,'units','characters','fontname','FixedWidth','fontsize',12);
value = gds.extravec(2);
rad  = uicontrol(handles.figuur,'Style','radiobutton','String','eps0','Tag','eps0par','Callback','set_option','Max',1,'Value',value,'BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(edit,'Callback','Hom_Callback');
user.pos = pose; set(edit,'Position',pose,'UserData',user);
    pos = [2 j 18 1.8]; user.pos=pos;
    set(rad,'Position',pos,'UserData',user); guidata(handles.figuur,handles);

j = j-2;
post = [2 j 18 1.8]; user.num = 0; user.pos = post; pose = [20 j 25 1.8];
tag  = strcat('edit',num2str(3));
edit = uicontrol(handles.figuur,'Style','edit','HorizontalAlignment','left','String',gds.eps1,'Tag',tag,'Backgroundcolor',color,'units','characters','fontname','FixedWidth','fontsize',12);
value = gds.extravec(3);
rad  = uicontrol(handles.figuur,'Style','radiobutton','String','eps1','Tag','eps1par','Callback','set_option','Max',1,'Value',value,'BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(edit,'Callback','Hom_Callback');
user.pos = pose; set(edit,'Position',pose,'UserData',user);
    pos = [2 j 18 1.8]; user.pos=pos;
    set(rad,'Position',pos,'UserData',user); guidata(handles.figuur,handles);

%--------------------------------------------------------------------------
function j = numeric_hetpars(handles,j)
color = [1 1 1];
j = j-2;
stat = uicontrol(handles.numeric_fig,'Style','text','String','Heteroclinic parameters','BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
pos = [5 j 38 1.8];user.num=0;user.pos=pos;
set(stat,'Position',pos,'UserData',user);
j = j-2;
post = [2 j 18 1.8];user.pos=post;
pose = [20 j 21 1.8];
stat1 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','T','Tag','text_T','BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
set(stat1,'Position',post,'UserData',user);user.pos=pose;
stat2 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','','Tag','T','Backgroundcolor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
set(stat2,'Position',pose,'UserData',user);
j = j-2;
post = [2 j 18 1.8];user.pos=post;
pose = [20 j 21 1.8];
stat1 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','eps0','Tag','text_eps0','BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
set(stat1,'Position',post,'UserData',user);user.pos=pose;
stat2 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','','Tag','eps0','Backgroundcolor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
set(stat2,'Position',pose,'UserData',user);
j = j-2;
post = [2 j 18 1.8];user.pos=post;
pose = [20 j 21 1.8];
stat1 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','eps1','Tag','text_eps1','BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
set(stat1,'Position',post,'UserData',user);user.pos=pose;
stat2 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','','Tag','eps1','Backgroundcolor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
set(stat2,'Position',pose,'UserData',user);
guidata(handles.numeric_fig,handles);

%--------------------------------------------------------------------------
function j = start_eigenvalues(handles,j)
global gds MC 
color = [1 1 1];
j = j-2;
pos = [12 j 18 1.8];user.num=0;user.pos=pos;
stat = uicontrol(handles.numeric_fig,'Style','text','String','Eigenvalues','Tag','eigenvalues','BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
set(stat,'Position',pos,'UserData',user);
jim = j -4*gds.dim;
for k = 1:2*gds.dim
    j = j-2;    
    tag1 = sprintf('Re[%d]',k);
    tag11 = sprintf('Re_%d',k);
    tag2 = sprintf('Im[%d]',k);
    tag22 = sprintf('Im_%d',k);
    post = [1 j 18 1.8];user.pos=post;
    pose = [19 j 27 1.8];             
    stat1 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String',tag1,'BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
    stat2 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','','Tag',tag11,'BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
    set(stat1,'Position',post,'UserData',user);user.pos=pose;
    set(stat2,'Position',pose,'UserData',user);
    jim=jim-2;
    post = [1 jim 18 1.8];user.pos=post;
    pose = [19 jim 27 1.8];             
    stat1 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String',tag2,'BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
    stat2 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','','Tag',tag22,'BackGroundColor',color,'fontname','FixedWidth','fontsize',12,'Units','characters');
    set(stat1,'Position',post,'UserData',user);user.pos=pose;
    set(stat2,'Position',pose,'UserData',user);

end
j=jim;