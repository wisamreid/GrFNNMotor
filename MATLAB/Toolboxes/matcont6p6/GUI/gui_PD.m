function  PD
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
global gds path_sys MC
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
        gds.point = 'PD ';
    end
    gds.type = 'PD ';
    set(0,'ShowHiddenHandles','on');        
    h = findobj('Type','uimenu','Tag','forward'); set(h,'enable','on');
    h = findobj('Type','uimenu','Tag','backward'); set(h,'enable','on');
    set(0,'ShowHiddenHandles','off');  
    matcont('make_curve');
    load_matcont;
    file = strcat(gds.system,'.mat');
    file = fullfile(path_sys,file);
    save(file,'gds');
    continuer; starter;
    if  ~isempty(MC.numeric_fig),numeric;end
    set(0,'ShowHiddenHandles','on');        
    h = findobj('Type','uimenu','Tag','window');set(h,'enable','on');
    h = findobj('Type','uimenu','Tag','extend'); set(h,'enable','off');
    set(0,'ShowHiddenHandles','off');        
%-------------------------------------------------------------------------    
function starter1(handles)
global gds  MC
    ndim = size(gds.parameters,1);
    color = [1 1 1];
    s = ndim*2+20+2*7+2;
    if isfield(gds,'userfunction')
        s = s+size(gds.userfunction,2)*2+2;
    end
    slider_step(1) = 2/s;
    slider_step(2) = 2/s;
    gds.options.ActiveUParams=[];
    gds.options.ActiveSParams=[];
    gds.options.ActiveSParam=[];    
    set(handles.slider1,'sliderstep',slider_step,'max',s,'min',0,'Value',s);
    j = s-2;
    stat = uicontrol(handles.figuur,'Style','text','String','Initial Point','Tag','Initial_Point','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
    pos = [10 j 35 1.80]; user.num = 0; user.pos = pos;
    set(stat,'Position',pos,'UserData',user);
    pos1 = starter('start_parameters',handles,j);
    pos1 = starter('start_period',handles,pos1);
    pos1 = starter('start_jacobian',handles,pos1);
    pos1 = starter('start_discret',handles,pos1);
    gds.options.IgnoreSingularity = 1:4;
    pos1 = start_monitor(handles,pos1);
    if isfield(gds,'userfunction'),pos1 = starter('start_userfunctions',handles,pos1);end
    starter('start_multipliers',handles,pos1);
    starter('in_start',handles);
%-------------------------------------------------------------------------    
function load_layout(handles)
global gds 
    if size(gds.numeric.PD,1)<7
         gds.numeric.PD{6,1} = 'user functions';
         gds.numeric.PD{6,2} = 0;
         gds.numeric.PD{7,1} = 'npoints';
         gds.numeric.PD{7,2} = 0;
     end
     for i = 1:7
         if gds.numeric.PD{i,2}==1
             gds.numeric.PD{i,1} = upper(gds.numeric.PD{i,1});
         elseif gds.numeric.PD{i,2}==0
             gds.numeric.PD{i,1} = lower(gds.numeric.PD{i,1});
         end
        string{i,1} = gds.numeric.PD{i,1};
    end    
    set(handles.layoutbox,'String',string);
%-------------------------------------------------------------------------    
function layoutbox(list,index_selected)
global gds 
    c = gds.numeric.PD{index_selected,2};    
    gds.numeric.PD{index_selected,2} = 1-c;    
%-------------------------------------------------------------------------    
function numeric1(handles)
global gds  MC
    ndim = size(gds.parameters,1);
     if size(gds.numeric.PD,1)<7
         gds.numeric.PD{6,1} = 'user functions';
         gds.numeric.PD{6,2} = 0;
         gds.numeric.PD{7,1} = 'npoints';
         gds.numeric.PD{7,2} = 0;
     end
    if isfield(gds,'userfunction')
        dimu = size(gds.userfunction,2);
    else dimu=0; 
    end
    s = 2;
    if gds.numeric.PD{1,2}==1
        s = s+2*ndim+2;
    end
    if gds.numeric.PD{2,2}==1
        s = s+4;
    end
    if gds.numeric.PD{3,2}==1
        s = s+8+2;
    end  
    if gds.numeric.PD{4,2}==1
        s = s+4*gds.dim+2;
    end  
    if gds.numeric.PD{5,2}==1
        s = s+2;
    end 
    if gds.numeric.PD{6,2}==1
        s = s+2*dimu+2;
    end
    if gds.numeric.PD{7,2}==1
        s = s+2;
    end    
    slider_step(1) = 2/s;
    slider_step(2) = 2/s;
    set(handles.slider1,'sliderstep',slider_step,'max',s,'min',0,'Value',s);
    j = s;    
    if gds.numeric.PD{1,2}==1
        j = numeric('start_parameters',handles,j);
    end
    if gds.numeric.PD{2,2}==1
        j = numeric('start_period',handles,j);
    end
    if gds.numeric.PD{3,2}==1
        j = start_testfunctions(handles,j);
    end
    if gds.numeric.PD{4,2}==1
        j = numeric('start_multipliers',handles,j);
    end  
    if gds.numeric.PD{5,2}==1
        j = numeric('start_stepsize',handles,j);
    end 
    if gds.numeric.PD{6,2}==1
        j = numeric('start_userfunctions',handles,j);
    end
    if gds.numeric.PD{7,2}==1
        j = numeric('start_npoints',handles,j);
    end

    numeric('in_numeric',handles);
%------------------------------------------------------------------------- 
 function [string,plotopts] = load_draw%plot-window
 global gds
    plotopts = {'Standard'; 'Maximum'; 'Minimum'; 'Norm'};
    string = {'Coordinates';'Parameters';'Period'};
    if gds.options.Multipliers==1
        string{4} = 'Multiplier';
    end   
    if ~isempty(gds.options.Userfunctions) && gds.options.Userfunctions ~= 0
        string{end+1} = 'Userfunction';
    end
%-------------------------------------------------------------------------    
function string = load_attributes
global gds 
    d=1;
    for k = 1:gds.dim
        if k < 10
            string{d,1} = sprintf('Mod[ %d ]',k);
            string{d+1,1} = sprintf('Arg[ %d ]',k);
        else
            string{d,1} = sprintf('Mod[ %d]',k);
            string{d+1,1} = sprintf('Arg[ %d]',k);
        end
        d = d+2;
    end
    string{d,1} = sprintf('Mod[All]');
    string{d+1,1} = sprintf('Arg[All]');
%-------------------------------------------------------------------------
function e = make_ready(pl,x)
global gds 
    ndim = size(gds.parameters,1);
    len  = size(pl,2);
    leng = size(x,1);
    e(1,len) = 0;
    for j = 1:len
        switch pl(j).type
        case 'Coordinates'
            if ~isfield(pl(j),'plotopt') || isempty(pl(j).plotopt)
                pl(j).plotopt = 'Standard';
            end
            switch pl(j).plotopt
                case 'Standard'
                    e(1,j) = pl(j).val;
                case 'Maximum'
                    e(1,j) = -1000*ndim - pl(j).val;
                case 'Minimum'
                    e(1,j) = -2000*ndim - pl(j).val;
                case 'Norm'
                    e(1,j) = -3000*ndim - pl(j).val;
            end
        case 'Parameters'
            if pl(j).val==gds.options.ActiveParams(1,1)
                e(1,j) = leng-1;
            elseif pl(j).val==gds.options.ActiveParams(1,2)
                e(1,j) = leng;
            else
                e(1,j) = -pl(j).val;
            end
        case 'Multiplier'
            if gds.options.Multipliers==0
                errordlg('It isn''t possible to plot the multipliers(they are not being computed)!');
            else
                e(1,j) = -ndim-pl(j).val;
            end
        case 'Period'
            e(1,j) = leng-2;
        case 'Userfunction'
            e(1,j) = -100000*ndim - pl(j).val;
        otherwise
            for k = 1:len
                e(1,k)=inf;
                return
            end
        end 
    end
%-------------------------------------------------------------------------    
function [label,lab] = label1(plot2)
global gds path_sys cds MC
    ndim = size(gds.parameters,1);
    len = size(plot2,2); dimplot = 0;
    label{len} = '';
    lab{len} = '';
    for k = 1:len
        label{k}='empty';
        lab{k}='empty';
        if plot2(1,k)<0 && plot2(1,k)>=-ndim
            pl2  = -plot2(1,k);
            label{k}  = sprintf('gds.parameters{%d,2}*ones(lds.tps,p)',pl2);  
            lab{k}  = sprintf('gds.parameters{%d,2}',pl2);      
        end     
        if plot2(1,k)==0
            label{k}  = 'gds.time{1,2}*ones(lds.tps,p)';
            lab{k} = 'gds.time{1,2}';
        end
        if plot2(1,k) < -ndim && plot2(1,k)>-1000*ndim %case multiplier
            if mod(plot2(1,k)+ndim,2)==0 % case Arg
                if (-plot2(1,k)-ndim)/2 > gds.dim % case All                    
                    label{k}  = sprintf('sort(angle(fout(lds.ntst+1+%d:lds.ntst+1+%d,i)))',(-plot2(1,k)-ndim)/2-gds.dim,(-plot2(1,k)-ndim)/2-1);
                    lab{k}= {};
                    for ii=1:gds.dim
                        lab{k,ii} = sprintf('angle(fout(lds.ntst+1+%d,i))',(-plot2(1,k)-ndim)/2-gds.dim+ii-1);
                    end
                else
                    label{k}  = sprintf('angle(fout(lds.ntst+1+%d,i))',(-plot2(1,k)-ndim)/2);
                    lab{k} = sprintf('angle(fout(lds.ntst+1+%d,i))',(-plot2(1,k)-ndim)/2);
                end
            else % case Mod
                if (-plot2(1,k)-ndim+1)/2 > gds.dim % case All                  
                    label{k}  = sprintf('sort(abs(fout(lds.ntst+1+%d:lds.ntst+1+%d,i)))',(-plot2(1,k)-ndim+1)/2-gds.dim,(-plot2(1,k)-ndim+1)/2-1);
                    lab{k}= {};
                    for ii=1:gds.dim
                        lab{k,ii} = sprintf('abs(fout(lds.ntst+1+%d,i))',(-plot2(1,k)-ndim+1)/2-gds.dim+ii-1);
                    end
                else
                    label{k} = sprintf('abs(fout(lds.ntst+1+%d,i))',(-plot2(1,k)-ndim+1)/2);
                    lab{k} = sprintf('abs(fout(lds.ntst+1+%d,i))',(-plot2(1,k)-ndim+1)/2);
                end
            end
        elseif plot2(1,k) <= -100000*ndim % case Userfunction
                label{k}  = sprintf('hout(2+%d,i)',-plot2(1,k)-100000*ndim); 
                lab{k}  = sprintf('hout(2+%d,i)',-plot2(1,k)-100000*ndim); 
        elseif plot2(1,k)<=-3000*ndim % case Norm
            label{k}  = sprintf('norm(xout((0:lds.tps-1)*lds.nphase+%d,i))',-plot2(1,k)-3000*ndim);  
        elseif plot2(1,k)<=-2000*ndim % case Minimum
            label{k}  = sprintf('min(xout((0:lds.tps-1)*lds.nphase+%d,i))',-plot2(1,k)-2000*ndim);    
        elseif plot2(1,k)<=-1000*ndim % case Maximum
            label{k}  = sprintf('max(xout((0:lds.tps-1)*lds.nphase+%d,i))',-plot2(1,k)-1000*ndim);   
        end
        if strcmp(label{k},'empty')==1
            dimplot = dimplot+1;
            if plot2(1,k) <= gds.dim %case coordinate
%                 dimplot = dimplot+1;
                label{k}  = sprintf('xout((0:lds.tps-1)*lds.nphase+%d,i)',plot2(1,k));                
            else %case active parameter
                label{k}  = sprintf('ones(lds.tps,1)*xout(%d,i)',plot2(1,k)); 
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
global gds  MC lds
    ndim = size(gds.parameters,1);
    num = 1;
    if isfield(MC.numeric_handles,'param1')
        label{ndim} = '';
        for d = 1:ndim
            label{num} = sprintf('set(MC.numeric_handles.param%d,''String'',num2str(parameters(%d),''%%0.8g''))',d,d);
            num = num +1;
        end            
    end
    if isfield(MC.numeric_handles,'period')
        label{num}='set(MC.numeric_handles.period,''String'',num2str(xout(end-2,i(end)),''%0.8g''))';
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
    if numsing(1)~=0 && isfield(MC.numeric_handles,'R2')   
        label{num}=sprintf('set(MC.numeric_handles.R2,''String'',num2str(hout(%d,i(end)),''%%.8g''))',numsing(1)+dimu);
        num = num + 1;
    end
    if numsing(2)~=0 && isfield(MC.numeric_handles,'FPD') 
        label{num}=sprintf('set(MC.numeric_handles.FPD,''String'',num2str(hout(%d,i(end)),''%%.8g''))',numsing(2)+dimu);
        num = num+1;
    end
    if numsing(3)~=0 && isfield(MC.numeric_handles,'GPD') 
        label{num}=sprintf('set(MC.numeric_handles.GPD,''String'',num2str(hout(%d,i(end)),''%%.8g''))',numsing(3)+dimu);
        num = num+1;
    end
    if numsing(4)~=0 && isfield(MC.numeric_handles,'PDNS') 
        label{num}=sprintf('set(MC.numeric_handles.PDNS,''String'',num2str(hout(%d,i(end)),''%%.8g''))',numsing(4)+dimu);
        num = num+1;
    end
    if (gds.options.Multipliers==1)&& isfield(MC.numeric_handles,'Mod_1')%multipliers numeric window
        for k = 1:gds.dim
            % Compute index of multipliers in fout
            indexmulti = length(lds.msh);
            label{num}=sprintf('set(MC.numeric_handles.Mod_%d,''String'',abs(fout(%d,i(end))))',k,indexmulti+k);
            label{num+1}=sprintf('set(MC.numeric_handles.Arg_%d,''String'',(180/pi)*angle(fout(%d,i(end))))',k,indexmulti+k);
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
global gds lds MC
    i=1;
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
        alll = 0;
        plot2 = varargin{2};
        if (plot2==[inf inf]),return;end
        hold on; d = axis;
        skew = 0.01*[d(2)-d(1) d(4)-d(3)]; 
        watchon;   axis(d);
        plo2{1}  = 'empt'; plo2{2}  = 'empt'; 
        plo2s{1} = 'empt'; plo2s{2} = 'empt'; dimplot = 0;
        for k = 1:2
            if plot2(1,k)<0 && plot2(1,k)>=-ndim
                pl2 = -plot2(1,k);
                plo2{k} = gds.parameters{pl2,2}*ones(lds.tps,p) ; 
            end   
            if (plot2(1,k)==0)
                plo2{k} = gds.time{1,2}*ones(lds.tps,p);
            end
            if plot2(1,k) < -ndim && plot2(1,k)>-1000*ndim % case multiplier
                if mod(plot2(1,k)+ndim,2)==1 % case Mod
                    if (-plot2(1,k)-ndim+1)/2 <= gds.dim
                        plo2{k} = abs(f(lds.ntst+1+(-plot2(1,k)-ndim+1)/2,:));
                        alll = 0;
                    else % case All
                        alll = 1;
                        tmp = [];
                        for iii = 1:gds.dim
                            abs(f(lds.ntst+1+(-plot2(1,k)-ndim+1)/2-gds.dim+iii-1,:));
                            tmp(:,end+1) = abs(f(lds.ntst+1+(-plot2(1,k)-ndim+1)/2-gds.dim+iii-1,:))';
                        end
                        plo2{k} = tmp;
                    end
                else % case Arg
                    if (-plot2(1,k)-ndim+1)/2 <= gds.dim
                        alll = 0;
                        plo2{k} = angle(f(lds.ntst+1+(-plot2(1,k)-ndim)/2,:));
                    else % case All
                        alll = 1;
                        tmp = [];
%                     label{k}  = sprintf('angle(fout(lds.ntst+1+%d,i))',(-plot2(1,k)-ndim)/2);
                        for iii = 1:gds.dim
                            tmp(:,end+1) = angle(f(lds.ntst+1+(-plot2(1,k)-ndim)/2-gds.dim+iii-1,:))';
                        end
                        plo2{k} = tmp;
                    end
                end 
        elseif plot2(1,k) <= -100000*ndim % case Userfunction
                plo2{k}  = ones(lds.tps,1)*h(2-plot2(1,k)-100000*ndim,:);
            
            elseif plot2(1,k)<=-3000*ndim % case Norm coords
                      plo2{k} = norm(x((0:lds.tps-1)*lds.nphase+(-plot2(1,k)-3000*ndim),:));
                  elseif plot2(1,k)<=-2000*ndim % case Min coords
                      plo2{k} = min(x((0:lds.tps-1)*lds.nphase+(-plot2(1,k)-2000*ndim),:));
                  elseif plot2(1,k)<=-1000*ndim % case Max coords
                      plo2{k} = max(x((0:lds.tps-1)*lds.nphase+(-plot2(1,k)-1000*ndim),:));
            end
            if strcmp(plo2{k},'empt')==1
                if plot2(1,k) <= gds.dim 
                    dimplot = dimplot+1;
                    plo2{k} = x((0:lds.tps-1)*lds.nphase+plot2(1,k),:);
                else
                    plo2{k} = ones(lds.tps,1)*x(plot2(1,k),:);
                end
            end            
        end
        if dimplot == 0 
            plo2{1} = plo2{1}(1,:);
            if alll == 0
                plo2{2} = plo2{2}(1,:);
            end
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
                        line(plo2{1}(:,ii),tmp(ii,:),'LineStyle','-','Color','b');
                    end
                else    
                    line(plo2{1},plo2{2},'LineStyle','-','Color','b');
                end
            end
            if ~strcmp(plo2s{1},'empt')&&~strcmp(plo2s{2},'empt')
                tmp = plo2s{2}';
                for ii = 1:size(tmp,1)
                    line(plo2s{1}(:,ii),tmp(ii,:),'Linestyle','-','Marker','.','MarkerEdgeColor','r');
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
                plo3{k} = gds.parameters{pl2,2}*ones(lds.tps,p); 
            elseif plo(1,k) <= -100000*ndim % case Userfunction
                plo3{k}  = ones(lds.tps,1)*h(2-plo(1,k)-100000*ndim,:);
            
            elseif plo(1,k)<=-3000*ndim % case Norm coords
                      plo3{k} = norm(x((0:lds.tps-1)*lds.nphase+(-plo(1,k)-3000*ndim),:));
                  elseif plo(1,k)<=-2000*ndim % case Min coords
                      plo3{k} = min(x((0:lds.tps-1)*lds.nphase+(-plo(1,k)-2000*ndim),:));
                  elseif plo(1,k)<=-1000*ndim % case Max coords
                      plo3{k} = max(x((0:lds.tps-1)*lds.nphase+(-plo(1,k)-1000*ndim),:));
            elseif plo(1,k) < -ndim % Case multipliers
                if mod(plo(1,k)+ndim,2)==1
                    plo3{k} = ones(lds.tps,1)*abs(f((-plo(1,k)-ndim+1)/2,:));
                else %Mod
                    plo3{k} = ones(lds.tps,1)*angle(f((-plo(1,k)-ndim)/2,:));
                end
            end     
            if (plo(1,k)==0) %case time
                plo3{k} = gds.time{1,2}*ones(lds.tps,p);
            end
            if strcmp(plo3{k},'empt')==1
                if plo(1,k)<=gds.dim       %case coordinate
                    dimplot = dimplot+1;
                    plo3{k}=x((0:lds.tps-1)*lds.nphase+plo(1,k),:);
                else %case active parameter
                    plo3{k}=x(plo(1,k)*ones(lds.tps,1),:);  
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
                line(plo3s{1},plo3s{2},plo3s{3},'linestyle','-','Marker','.','MarkerEdgecolor','r');
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
        xout=x;hout=h;fout=f;
        for d=1:ndim
            parameters(d) = gds.parameters{d,2};
        end   
        j = gds.options.ActiveParams(1,1);
        parameters(j) = xout((gds.dim)+1,s.index);
        dat = get(MC.numeric_fig,'Userdata');
        i = s.index;
        npoints=s.index;
        cds.h='?';
        for k=1:size(dat.label,2)
            eval(dat.label{k});  
        end 
    end  
%-------------------------------------------------------------------------
function output(numsing,xout,s, hout,fout,i)
global gds path_sys cds lds MC
    ndim = size(gds.parameters,1);
    x1 = xout(:,i);npoints=i(end);
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
                        line(s1(1,:),s2(iii,:),'Color','b','Linestyle','-'); 
                  end             
                else
                    line(s1,s2,'Color','b','Linestyle','-'); 
                end
%                 line(s1(1,:),s2(1,:),'Color','b','Linestyle','-'); 
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
                line(s1(1,:),s2(1,:),s3(1,:),'linestyle','-','color','b');                  
            end
            drawnow;
        end         
    end
    if ~isempty(MC.numeric_fig) %numeric window is open
        for d = 1:ndim
            parameters(d) = gds.parameters{d,2};
        end   
        for k = 1:2           
            j = gds.options.ActiveParams(1,k);         
            parameters(j) = xout(end-2+k,i(end));                
        end        
        dat = get(MC.numeric_fig,'Userdata');
        for k=1:size(dat.label,2)
            eval(dat.label{k});  
        end
    end
%-------------------------------------------------------------------------    
function load_point(index,x,string,file,varargin)
global gds lds 
    load(file);
    for i = 1:gds.dim
        gds.coordinates{i,2} = x(i,index);
    end
    gds.discretization.ntst = lds.ntst;
    gds.discretization.ncol = lds.ncol;
    ndim = size(lds.P0,1);
    for i = 1:ndim
        gds.parameters{i,2} = lds.P0(i,1);
    end    
    len = length(lds.ActiveParams);  
    if len>0
        for i = 1:len
            j = lds.ActiveParams(1,i);
            gds.parameters{j,2} = x(end-i+1,index);
        end
    end
    gds.period = x(end-len,index);
    if strcmp(string,'save')
        num = varargin{1};
        save(file,'x','v','s','h','f','num','cds','lds','ctype','point');
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
global gds path_sys cds 
     if ~isempty(varargin)
        file = fullfile(path_sys,gds.system,gds.diagram,strcat(gds.curve.new,'.mat'));
        if exist(file,'file' )
            load(file);
        else
            errordlg('It is not possible to extend the current curve!','Error extend PD');
            return 
        end  
        [x,v,s,h,f] = cont(x,v,s,h,f,cds);
     else
        if (size(gds.options.ActiveParams,2)~=2)
            errordlg('You have to select exactly two free parameters');
            return;
        end    
        if ~isempty(gds.curve.old)
            file =  strcat(gds.curve.old,'.mat');
            file = fullfile(path_sys,gds.system,gds.diagram,file);
            if exist(file,'file' )
                load(file);        
            else
                errordlg('not enough information available','Error start cont PD');
                return
            end
        else
            errordlg('not enough information available','Error start cont PD');
            return
        end
        systemhandle = str2func(gds.system);
        switch gds.point
        case 'GPD '
            [x0,v0] = init_GPD_PD(systemhandle,x,s(num),gds.options.ActiveParams,gds.discretization.ntst,gds.discretization.ncol);    
        case 'LPPD '
            [x0,v0] = init_LPPD_PD(systemhandle,x,s(num),gds.options.ActiveParams,gds.discretization.ntst,gds.discretization.ncol);    
        case 'PDNS '
            [x0,v0] = init_PDNS_PD(systemhandle,x,s(num),gds.options.ActiveParams,gds.discretization.ntst,gds.discretization.ncol);    
        case 'R2 '
            [x0,v0] = init_R2_PD(systemhandle,x,s(num),gds.options.ActiveParams,gds.discretization.ntst,gds.discretization.ncol);    
        otherwise
            [x0,v0] = init_PD_PD(systemhandle,x,s(num),gds.options.ActiveParams,gds.discretization.ntst,gds.discretization.ncol);    
        end
        if isempty(x0)
            return
        end       
        [x,v,s,h,f] = cont(@perioddoubling,x0,v0,gds.options); 
    end
    file = fullfile(path_sys,gds.system,gds.diagram,gds.curve.new);
    ctype  = gds.type; point = gds.point;
    status = mkdir(path_sys,gds.system);
    dir = fullfile(path_sys,gds.system);
    status = mkdir(dir,gds.diagram);
    if ~exist('num'),num = 1;end
    if ~isempty(x)
        save(file,'x','v','s','h','f','cds','lds','ctype','point','num');          
    end


%---------------------------------------------------------------------------
function j = start_monitor(handles,j)
global gds
color = [1 1 1];
j = j-2;
pos = [10 j 35 1.80]; user.num = 0; user.pos = pos;
stat = uicontrol(handles.figuur,'Style','text','String','Monitor Singularities','Tag','singularities','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',pos,'UserData',user);
%BP=1;H=1;LP=3
for k = 1:4    
    tag1 = strcat('text',num2str(k));           
    switch k
    case 1
        if gds.dim < 3
            continue
        end
        string = '1:2 resonance';
    case 2
        if gds.dim < 3
            continue
        end
        string = 'Fold-flip';
    case 3
        if gds.dim < 2
            continue
        end       
        string = 'Generalized PD';
    case 4       
        if (gds.dim <4)
            continue;
        end
        string = 'Flip-NS';            
    otherwise
        errordlg('Is not possible');
    end
    if ~isempty(gds.options.IgnoreSingularity)&find(gds.options.IgnoreSingularity==k),valu=2;else valu=1;end
    j = j-2; user.num = num2str(k);
    post = [2 j 20 1.80];user.pos = post;
    pose = [22 j 23 1.80];
    stat = uicontrol(handles.figuur,'Style','text','HorizontalAlignment','left','String',string,'Tag',tag1,'BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
    edit = uicontrol(handles.figuur,'Style','popupmenu','String',{'yes','no'},'Tag','IgnoreSingularity','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12,'Value',valu);
    set(stat,'Position',post,'UserData',user); user.pos = pose;
    set(edit,'Position',pose,'Callback','set_option','UserData',user);
end
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
%R2;LPPD;GPD;PDNS
j = j-2;
post = [2 j 20 1.80];user.pos=post;
pose = [22 j 22 1.80];
%'R2';
if gds.dim > 2    
    stat1 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','R2','BackGroundColor',color,'Units','characters','fontname','FixedWidth','fontsize',12);
    stat2 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','','Tag','R2','BackGroundColor',color,'Units','characters','fontname','FixedWidth','fontsize',12);
    set(stat1,'Position',post,'UserData',user);user.pos=pose;
    set(stat2,'Position',pose,'UserData',user);
    j = j-2;
    post = [2 j 20 1.80];user.pos=post;
    pose = [22 j 22 1.80];       
    %LPPD
    stat1 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','LPPD','BackGroundColor',color,'Units','characters','fontname','FixedWidth','fontsize',12);
    stat2 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','','Tag','LPPD','BackGroundColor',color,'Units','characters','fontname','FixedWidth','fontsize',12);
    set(stat1,'Position',post,'UserData',user);user.pos=pose;
    set(stat2,'Position',pose,'UserData',user);
    j = j-2;
    post = [2 j 20 1.80];user.pos=post;
    pose = [22 j 22 1.80];       
end
%GPD
if gds.dim > 1
    stat1 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','GPD','BackGroundColor',color,'Units','characters','fontname','FixedWidth','fontsize',12);
    stat2 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','','Tag','GPD','BackGroundColor',color,'Units','characters','fontname','FixedWidth','fontsize',12);
    set(stat1,'Position',post,'UserData',user);user.pos=pose;
    set(stat2,'Position',pose,'UserData',user);
    j = j-2;
    post = [2 j 20 1.80];user.pos=post;
    pose = [22 j 22 1.80];       
end
%PDNS
if gds.dim > 3    
    stat1 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','PDNS','BackGroundColor',color,'Units','characters','fontname','FixedWidth','fontsize',12);
    stat2 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','','Tag','PDNS','BackGroundColor',color,'Units','characters','fontname','FixedWidth','fontsize',12);
    set(stat1,'Position',post,'UserData',user);user.pos=pose;
    set(stat2,'Position',pose,'UserData',user);
end
guidata(handles.numeric_fig,handles);

