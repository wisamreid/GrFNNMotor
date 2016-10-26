function gui_EP
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
    
%------------------------------------------------------
function point(tag)
global gds path_sys MC
    %tagexample:EP_EP(it is only possible to start a
                      % EP-curve continuation from a EP)
    [point,str] = strtok(tag,'_');%point=EP;str=EP
    set(0,'ShowHiddenHandles','on');
    h = findobj('Type','uimenu','Tag','window');
    set(h,'enable','on');
    h = findobj('Type','uimenu','Tag','curvetype');
    set(0,'ShowHiddenHandles','off');
    list = get(h,'children');%find all the type of curves
    tag  = get(list,'Tag'); 
    str  = strcat(str,'_');
    for i = 1:length(list)
        if isempty(findstr(str,strcat(tag{i},'_')))
            %makes other continuations of curves not possible
            set(list(i),'enable','off');
        else
            %EP-curve is possible
            set(list(i),'enable','on');
        end
    end     
    %gds.point='EP '
    gds.point = sprintf('%s%c',point,char(32)); 
    %callback of 'Type->curve->equilibrium'
    if strmatch('__',str,'exact')
        set(MC.mainwindow.compute,'enable','off');
        gds.type='';
        gds.curve.new='';
    else curvetype;
    end
%------------------------------------------------------   
function curvetype
global gds path_sys MC
   if isempty(gds.point)
        gds.point = 'EP ';%the default
    end
    gds.type = 'EP ';
    set(0,'ShowHiddenHandles','on');        
    h = findobj('Type','uimenu','Tag','forward'); set(h,'enable','on');
    h = findobj('Type','uimenu','Tag','backward'); set(h,'enable','on');
    h = findobj('Type','uimenu','Tag','extend'); set(h,'enable','on');
    set(0,'ShowHiddenHandles','off');  
    %makes curvename: EP_EP(number)
    matcont('make_curve');
    %load the new settings in the main matcont-window
    load_matcont;
    %save the new settings in the system file
    file = strcat(gds.system,'.mat');
    file = fullfile(path_sys,file);
    save(file,'gds');
    %load the appropriate continuer and starter window (see starter)
    continuer;starter;
    if  ~isempty(MC.numeric_fig),numeric;end
    %make it possible to open the plot/numeric window
    set(0,'ShowHiddenHandles','on');        
    h = findobj('Type','uimenu','Tag','window');
    set(h,'enable','on');
    set(0,'ShowHiddenHandles','off');  
%------------------------------------------------------
function starter1(handles)
global gds MC
    ndim=size(gds.parameters,1);
    %default is white
    color = [1 1 1];
    %the size of the window 2/unit
    s = 20+2*gds.dim+ndim*2;
    if isfield(gds,'userfunction')
        s=s+size(gds.userfunction,2)*2+2;
    end
    gds.options.ActiveUParams=[];
    gds.options.ActiveSParams=[];
    gds.options.ActiveSParam=[];    
    if strcmp(gds.point,'BP '),s = s+4;end
    slider_step(1) = 2/s;
    slider_step(2) = 2/s;
    set(handles.slider1,'sliderstep',slider_step,'max',s,'min',0,'Value',s);
    j = s-2;
    stat = uicontrol(handles.figuur,'Style','text','String','Initial Point','Tag','Initial_Point','BackGroundColor',color,'Units','characters','fontname','FixedWidth','fontsize',12);
    pos = [10 j 35 1.80];user.pos=pos;user.num=0;
    set(stat,'Position',pos,'UserData',user);     
    %make place for coordinates(2+2*gds.dim)
    pos1 = starter('start_coordinates',handles,j);
    %make place for parameters(ndim*2)
    pos1 = starter('start_parameters',handles,pos1);
    %make place for jacobian-settings(2+2)
    pos1 = starter('start_jacobian',handles,pos1);
    if isfield(gds,'userfunction'),pos1 = starter('start_userfunctions',handles,pos1);end
    if strcmp(gds.point,'BP '),pos1=start_amplitude(handles,pos1);end
    %make place for testfunctions(2+3*2)
    pos1 = start_monitor(handles,pos1);
    %make place for eigenvalues(2*2)
    pos1 = starter('start_eigenvalues',handles,pos1);
    starter('in_start',handles);
%------------------------------------------------------    
function load_layout(handles)
global gds path_sys
    ndim=size(gds.parameters,1);
    if size(gds.numeric.EP,1)<7
        gds.numeric.EP{6,1} = 'user functions';
        gds.numeric.EP{6,2} = 0;
        gds.numeric.EP{7,1} = 'npoints';
        gds.numeric.EP{7,2} = 0;
    end
    for i =1:size(gds.numeric.EP,1)
        if gds.numeric.EP{i,2}==1
            gds.numeric.EP{i,1} = upper(gds.numeric.EP{i,1});
        elseif gds.numeric.EP{i,2}==0
            gds.numeric.EP{i,1} = lower(gds.numeric.EP{i,1});
        end
    end
    file = fullfile(path_sys,gds.system);
    save(file,'gds');
    string = cellstr(char(gds.numeric.EP{:,1}));
    set(handles.layoutbox,'String',string);
%-----------------------------------------------------
function layoutbox(list,index_selected)    
global gds
    b = gds.numeric.EP{index_selected,2}; 
    gds.numeric.EP{index_selected,2} = 1-b;
%-----------------------------------------------------
function numeric1(handles)
global gds 
    ndim=size(gds.parameters,1);
    if size(gds.numeric.EP,1)<7
        gds.numeric.EP{6,1} = 'user functions';
        gds.numeric.EP{6,2} = 0;
        gds.numeric.EP{7,1} = 'npoints';
        gds.numeric.EP{7,2} = 0;
    end
   if isfield(gds,'userfunction')
        dimu = size(gds.userfunction,2);
    else dimu=0; 
   end
    s=2;
    if gds.numeric.EP{6,2}==1
        s=s+2*dimu+2;
    end
    if gds.numeric.EP{1,2}==1
        s =s+2*gds.dim+2;
    end 
    if gds.numeric.EP{2,2}==1        
        s=s+2*ndim+2;
    end
    if gds.numeric.EP{3,2}==1
        s=s+3*2+2;
    end
    if gds.numeric.EP{4,2}==1
        s=4*gds.dim+2+s;
    end
    if gds.numeric.EP{5,2}==1
        s=s+3;
    end
    if gds.numeric.EP{7,2}==1
        s=s+2;
    end
    slider_step(1) = 2/s;
    slider_step(2) = 2/s;
    set(handles.slider1,'sliderstep',slider_step,'max',s,'min',0,'Value',s);
    j = s;
    if gds.numeric.EP{1,2}==1
        j = numeric('start_coordinates',handles,j);
    end 
    if gds.numeric.EP{2,2}==1        
        j = numeric('start_parameters',handles,j);
    end
    if gds.numeric.EP{3,2}==1
        j = start_testfunctions(handles,j);
    end
    if gds.numeric.EP{4,2}==1
        j = numeric('start_eigenvalues',handles,j);
    end
    if gds.numeric.EP{5,2}==1
        j = numeric('start_stepsize',handles,j);
    end
    if gds.numeric.EP{6,2}==1
        j = numeric('start_userfunctions',handles,j);
    end
    if gds.numeric.EP{7,2}==1
        j = numeric('start_npoints',handles,j);
    end
    numeric('in_numeric',handles);
%------------------------------------------------------
function [string, plotopts] = load_draw
global gds 
    plotopts = {};
    string ={'Coordinates';'Parameters'};
    if gds.options.Multipliers==1
        string{3} = 'Eigenvalues';
    end
    if ~isempty(gds.options.Userfunctions) && gds.options.Userfunctions ~= 0
        string{end+1} = 'Userfunction';
    end
%------------------------------------------------------
function string = load_attributes%plot-window
global gds 
if gds.dim < 10
    d = 1;
    for k = 1:gds.dim
        string{d,1} = sprintf(' %d ',k);
        d = d + 1;  
    end
    string{d,1} = 'All';
else
    d = 1;
    for k = 1:gds.dim
        if k<10
            string{d,1} = sprintf(' 0%d',k);
        else
            string{d,1} = sprintf(' %d',k);
        end
        d = d + 1;  
    end
    string{d,1} = 'All';
end
    
%------------------------------------------------------    
function e = make_ready(pl,x)
global gds
    ndim=size(gds.parameters,1);
    len = size(pl,2);
    e(1,len) = 0;
    for j = 1:len
        switch pl(j).type
        case 'Coordinates'
            e(1,j) = pl(j).val;
        case 'Parameters'
            p = ismember(pl(j).val,gds.options.ActiveParams);
            if p
                e(1,j) = 1+gds.dim;
            else
                e(1,j) = -pl(j).val;
            end
        case 'Time'
            e(1,j) = 0;
        case 'Eigenvalues'
            if gds.options.Eigenvalues==0
                errordlg('It isn''t possible to plot the multipliers (they are not being computed)!');
            else
                switch pl(j).eigsopt
                    case 'Real part'
                        e(1,j) = -1000*ndim-pl(j).val;
                    case 'Imag. part'
                        e(1,j) = -2000*ndim-pl(j).val;
                    case 'Norm'
                        e(1,j) = -3000*ndim-pl(j).val;
                end
                if pl(j).val > gds.dim
                    e(1,j) = e(1,j) -10000*ndim;
                end
            end
        case 'Userfunction'
            e(1,j) = -100000*ndim - pl(j).val;
        otherwise
            for k = 1:len
                e(1,k)=inf;
                return
            end
        end 
    end    
%------------------------------------------------------
function [label,labels] = label1(plot2)
global gds 
    ndim = size(gds.parameters,1);
    len = size(plot2,2);
    label{len} = '';
    labels{len} = '';
    for k = 1:len
        label{k}='empty';
        if plot2(1,k)<0 && plot2(1,k)>=-ndim
            pl2  = -plot2(1,k);
            label{k}  = sprintf('gds.parameters{%d,2}*ones(1,ps)',pl2);
            labels{k} = sprintf('gds.parameters{%d,2}',pl2); 
        end     
        if plot2(1,k)==0
            label{k}  ='gds.time{1,2}*ones(1,ps)';
            labels{k} ='gds.time{1,2}';
        end
        if plot2(1,k) < -ndim %case eigenvalues
            plotall = 0;
            if plot2(1,k) <= -100000*ndim % case Userfunction
                label{k}  = sprintf('hout(2+%d,i(1)-jj:i(end))',-plot2(1,k)-100000*ndim); 
                labels{k}  = sprintf('hout(2+%d,i(sind))',-plot2(1,k)-100000*ndim); 
            elseif plot2(1,k) < -10000*ndim
                plotall = 1;
                plot2(1,k) = plot2(1,k) + 10000*ndim;
            elseif plot2(1,k) <= -3000*ndim
                if ~plotall
                    label{k}  = sprintf('abs(fout(%d,i(1)-jj:i(end)))',-plot2(1,k)-3000*ndim);
                    labels{k} = sprintf('abs(fout(%d,i(sind)))',-plot2(1,k)-3000*ndim);
                else
                    label{k}  = sprintf('abs(fout(%d:%d,i(1)-jj:i(end)))',-plot2(1,k)-3000*ndim-gds.dim,-plot2(1,k)-3000*ndim-1);
                    labels{k}= [];
                    for ii=1:gds.dim
                        tmp = labels{k};
                        tmp(ii,:) = sprintf('abs(fout(%d,i(sind)))',-plot2(1,k)-3000*ndim-gds.dim+ii-1);
                        labels{k} = char(tmp);
                    end
                end
            elseif plot2(1,k) <= -2000*ndim %Imag
                if ~plotall
                    label{k}  = sprintf('imag(fout(%d,i(1)-jj:i(end)))',-plot2(1,k)-2000*ndim);
                    labels{k} = sprintf('imag(fout(%d,i(sind)))',-plot2(1,k)-2000*ndim);
                else
                    label{k}  = sprintf('imag(fout(%d:%d,i(1)-jj:i(end)))',-plot2(1,k)-2000*ndim-gds.dim,-plot2(1,k)-2000*ndim-1);
                    labels{k}= [];
                    for ii=1:gds.dim
                        tmp = labels{k};
                        tmp(ii,:) =  sprintf('imag(fout(%d,i(sind)))',-plot2(1,k)-2000*ndim-gds.dim+ii-1);
                        labels{k} = char(tmp);
                    end
                end
            else %Re
                if ~plotall
                    label{k} = sprintf('real(fout(%d,i(1)-jj:i(end)))',-plot2(1,k)-1000*ndim);
                    labels{k}= sprintf('real(fout(%d,i(sind)))',-plot2(1,k)-1000*ndim);
                else
                    label{k} = sprintf('real(fout(%d:%d,i(1)-jj:i(end)))',-plot2(1,k)-1000*ndim-gds.dim,-plot2(1,k)-1000*ndim-1);
                    labels{k}= [];
                    for ii=1:gds.dim
                        tmp = labels{k};
                        tmp(ii,:) =  sprintf('real(fout(%d,i(sind)))',-plot2(1,k)-1000*ndim-gds.dim+ii-1);
                        labels{k} = char(tmp);
                    end
                end
            end
        end
        if strcmp(label{k},'empty')==1
            label{k} = sprintf('xout(%d,i(1)-jj:i(end))',plot2(1,k));
            labels{k} = sprintf('xout(%d,i(sind))',plot2(1,k));
        end   
    end 
%------------------------------------------------------
function label = numeric_label(numsing) 
global gds MC
    ndim=size(gds.parameters,1);
    num = 1;
    label{gds.dim} = '';
    if isfield(MC.numeric_handles,'coord1')
        for d = 1:gds.dim
            label{num}= sprintf('set(MC.numeric_handles.coord%d,''String'',num2str(xout(%d,i(end)),''%%0.8g''))',d,d);    
            num = num+1;
        end        
    end
    if isfield(MC.numeric_handles,'param1')
        for d = 1:ndim
            label{num} = sprintf('set(MC.numeric_handles.param%d,''String'',num2str(parameters(%d),''%%0.8g''))',d,d);
            num = num +1;
        end            
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
    if numsing(1)~=0 && isfield(MC.numeric_handles,'BP')  
        label{num}=sprintf('set(MC.numeric_handles.BP,''String'',num2str(hout(%d,i(end)),''%%0.8g''))',numsing(1)+dimu);
        num = num + 1;
    end
    if numsing(2)~=0 && isfield(MC.numeric_handles,'H') 
        label{num}=sprintf('set(MC.numeric_handles.H,''String'',num2str(hout(%d,i(end)),''%%0.8g''))',numsing(2)+dimu);
        num = num+1;
    end
    if numsing(3)~=0 && isfield(MC.numeric_handles,'LP') 
        label{num} = sprintf('set(MC.numeric_handles.LP,''String'',num2str(hout(%d,i(end)),''%%0.8g''))',numsing(3)+dimu);
        num = num+1;
    end
    if (gds.options.Eigenvalues==1)&& isfield(MC.numeric_handles,'Re_1')%multipliers numeric window
        for k = 1:gds.dim
            label{num}=sprintf('set(MC.numeric_handles.Re_%d,''String'',real(fout(%d,i(end))))',k,k);
            label{num+1}=sprintf('set(MC.numeric_handles.Im_%d,''String'',imag(fout(%d,i(end))))',k,k);
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
%------------------------------------------------------
function draw(varargin)    
global gds path_sys cds eds MC
    ndim = size(gds.parameters,1);
    if strcmp(varargin{4},'select')
        x = varargin{3}; s = varargin{5};h = varargin{6}; f = varargin{7};
    else
        file = varargin{3};load(file);
        s(1)   = [];   s(end) = [];
    end
   switch varargin{1}
    case 2
        plot2 = varargin{2};
        if (plot2==[inf inf]),return;end
        hold on; d = axis;
        skew = 0.01*[d(2)-d(1) d(4)-d(3)];
        axis(d);plo2{1} = 'empt';plo2{2} = 'empt' ;    
        plo2s{1} = 'empt';plo2s{2} = 'empt';      
        p = size(x,2);
        k = size(cat(1,s.index),1);
        for j = 1:2
            if plot2(1,j)<0 && plot2(1,j)>=-ndim
                pl2 = -plot2(1,j);
                plo2{j}  = gds.parameters{pl2,2}*ones(1,p);
                plo2s{j} = gds.parameters{pl2,2}*ones(1,k);
            end  
            if (plot2(1,j)==0)%case time
                plo2{j}  = gds.time{1,2}*ones(1,p);
                plo2s{j} = gds.time{1,2}*ones(1,k);  
            end      
        if plot2(1,j) < -ndim %case eigenvalues
            
            plotall = 0;
            if plot2(1,j) <= -100000*ndim % case Userfunction
                plo2{j}  = h(2-plot2(1,j)-100000*ndim,:); 
                plo2s{j} =  h(2-plot2(1,j)-100000*ndim,cat(1,s.index)); 
            
            elseif plot2(1,j) < -10000*ndim
                plotall = 1;
                plot2(1,j) = plot2(1,j) + 10000*ndim;
            elseif plot2(1,j) <= -3000*ndim
                if ~plotall
                    plo2{j} = abs(f(-plot2(1,j)-3000*ndim,:));
                    plo2s{j}= abs(f(-plot2(1,j)-3000*ndim,cat(1,s.index)));
                else % case All
                    plo2{j} = abs(f(-plot2(1,j)-3000*ndim-gds.dim:-plot2(1,j)-3000*ndim-1,:))';
                    plo2s{j}= abs(f(-plot2(1,j)-3000*ndim-gds.dim:-plot2(1,j)-3000*ndim-1,cat(1,s.index)))';
                end
            elseif plot2(1,j) <= -2000*ndim %Imag
                if ~plotall
                    plo2{j} = imag(f(-plot2(1,j)-2000*ndim,:));
                    plo2s{j}= imag(f(-plot2(1,j)-2000*ndim,cat(1,s.index)));
                else % case All
                    plo2{j} = imag(f(-plot2(1,j)-2000*ndim-gds.dim:-plot2(1,j)-2000*ndim-1,:))';
                    plo2s{j}= imag(f(-plot2(1,j)-2000*ndim-gds.dim:-plot2(1,j)-2000*ndim-1,cat(1,s.index)))';
                end
            else %Re
                if ~plotall
                    plo2{j} = real(f(-plot2(1,j)-1000*ndim,:));
                    plo2s{j}= real(f(-plot2(1,j)-1000*ndim,cat(1,s.index)));
                else
                    plo2{j} = real(f(-plot2(1,j)-1000*ndim-gds.dim:-plot2(1,j)-1000*ndim-1,:))';
                    plo2s{j}= real(f(-plot2(1,j)-1000*ndim-gds.dim:-plot2(1,j)-1000*ndim-1,cat(1,s.index)))';
                end
            end
        end
            if strcmp(plo2{j},'empt')==1
                plo2{j}  = x(plot2(1,j),:);
                plo2s{j} = x(plot2(1,j),cat(1,s.index));
            end
        end  
        if strcmp(varargin{4},'redraw')      
            if ~strcmp(plo2{1},'empt')&&~strcmp(plo2{2},'empt')
              line(plo2{1},plo2{2},'LineStyle','-','Color','b');
            end
            if ~strcmp(plo2s{1},'empt')&&~strcmp(plo2s{2},'empt')
                line(plo2s{1},plo2s{2},'LineStyle','none' , 'Marker' , '*','Color','r');
            end     
            if ~isempty(s)       
                t1 = plo2s{1};
                t2 = plo2s{2};
                if length(t1) > length(t2)
                    t2 = t2*ones(size(t1));
                elseif length(t2) > length(t1)
                    t1 = t1*ones(size(t2));
                end                    
                text(t1+skew(1), t2+skew(2), cat(1,char(s.label)));            
            end
        else
            plot(plo2s{1},plo2s{2},'m+','MarkerSize',200);                
             hold off;
        end
    case 3 %3D
        plo = varargin{2};
        if (plo==[inf inf inf]),return;end
        hold on;d=axis;
        skew = 0.01*[d(2)-d(1) d(4)-d(3) d(6)-d(5)];
        p = size(x,2);k=size(cat(1,s.index),1);
        axis(d);plo3{1} = 'empt';plo3{2} = 'empt';plo3{3} = 'empt';
        plo3s{1} = 'empt'; plo3s{2} = 'empt';plo3s{3} = 'empt';
        for j = 1:3
            if plo(1,j)<0 && plo(1,j)>=-ndim
                pl2  = -plo(1,j);
                plo3{j}  = gds.parameters{pl2,2}*ones(1,p);
                plo3s{j} = gds.parameters{pl2,2}*ones(1,k);
            end       
            if (plo(1,j)==0)
                plo3{j}  = gds.time{1,2}*ones(1,p);
                plo3s{j} = gds.time{1,2}*ones(1,k);
            end
        if plo(1,j) < -ndim %case eigenvalues
            plotall = 0;
            if plo(1,j) <= -100000*ndim % case Userfunction
                plo3{j}  = h(2-plo(1,j)-100000*ndim,:); 
                plo3s{j} =  h(2-plo(1,j)-100000*ndim,cat(1,s.index)); 
            
            elseif plo(1,j) < -10000*ndim
                plotall = 1;
                plo(1,j) = plo(1,j) + 10000*ndim;
            elseif plo(1,j) <= -3000*ndim
                if ~plotall
                    plo3{j} = abs(f(-plo(1,j)-3000*ndim,:));
                    plo3s{j}= abs(f(-plo(1,j)-3000*ndim,cat(1,s.index)));
                else % case All
                    plo3{j} = abs(f(-plo(1,j)-3000*ndim-gds.dim:-plo(1,j)-3000*ndim-1,:))';
                    plo3s{j}= abs(f(-plo(1,j)-3000*ndim-gds.dim:-plo(1,j)-3000*ndim-1,cat(1,s.index)))';
                end
            elseif plo(1,j) <= -2000*ndim %Imag
                if ~plotall
                    plo3{j} = imag(f(-plo(1,j)-2000*ndim,:));
                    plo3s{j}= imag(f(-plo(1,j)-2000*ndim,cat(1,s.index)));
                else % case All
                    plo3{j} = imag(f(-plo(1,j)-2000*ndim-gds.dim:-plo(1,j)-2000*ndim-1,:))';
                    plo3s{j}= imag(f(-plo(1,j)-2000*ndim-gds.dim:-plo(1,j)-2000*ndim-1,cat(1,s.index)))';
                end
            else %Re
                if ~plotall
                    plo3{j} = real(f(-plo(1,j)-1000*ndim,:));
                    plo3s{j}= real(f(-plo(1,j)-1000*ndim,cat(1,s.index)));
                else
                    plo3{j} = real(f(-plo(1,j)-1000*ndim-gds.dim:-plo(1,j)-1000*ndim-1,:))';
                    plo3s{j}= real(f(-plo(1,j)-1000*ndim-gds.dim:-plo(1,j)-1000*ndim-1,cat(1,s.index)))';
                end
            end
        end
            if strcmp(plo3{j},'empt')==1
                plo3s{j} = x(plo(1,j),cat(1,s.index));
                plo3{j}  = x(plo(1,j),:);
            end
        end  
        if strcmp(varargin{4},'redraw')      
            if ~strcmp(plo3{1},'empt')&&~strcmp(plo3{2},'empt')&&~strcmp(plo3{3},'empt')
                line(plo3{1},plo3{2},plo3{3},'linestyle','-','Color','b');
            end
             if ~strcmp(plo3s{1},'empt')&&~strcmp(plo3s{2},'empt')&&~strcmp(plo3s{3},'empt')        
                t1 = plo3s{1};
                t2 = plo3s{2};
                t3 = plo3s{3};
                mx = max([length(t1),length(t2),length(t3)]);
                if mx > length(t1)
                    t1 = t1*ones(size(t1,1),mx);
                end
                if mx > length(t2)
                    t2 = t2*ones(size(t2,1),mx);
                end
                if mx > length(t3)
                    t3 = t3*ones(size(t3,1),mx);
                end              
                line(t1,t2,t3,'linestyle', 'none' , 'Marker' , '*','Color','r');
            end
            if ~isempty(s)             
                t1 = plo3s{1};
                t2 = plo3s{2};
                t3 = plo3s{3};
                mx = max([length(t1),length(t2),length(t3)]);
                if mx > length(t1)
                    t1 = t1*ones(size(t1,1),mx);
                end
                if mx > length(t2)
                    t2 = t2*ones(size(t2,1),mx);
                end
                if mx > length(t3)
                    t3 = t3*ones(size(t3,1),mx);
                end              
                text(t1+skew(1), t2+skew(2),t3+skew(3),cat(1,char(s.label)) );       
            end
        else            
            plot3(plo3s{1},plo3s{2},plo3s{3},'m+','MarkerSize',200);                
            hold off;
        end
    case 4 %3D
        xout=x;hout=h;fout=f;
        for d=1:ndim
            parameters(d) = gds.parameters{d,2};
        end   
        j = gds.options.ActiveParams(1,1);
        parameters(j) = xout((gds.dim)+1,s.index);
        dat = get(MC.numeric_fig,'Userdata');
        i=s.index;
        npoints=s.index;
        cds.h='?';
        for k=1:size(dat.label,2)
            eval(dat.label{k});  
        end 
   end  
%------------------------------------------------------    
function output(numsing,xout,s,hout,fout,i)
global gds MC cds eds
    ndim = size(gds.parameters,1);
    x1 = xout(:,i); npoints=i(end);
    if strcmp(s.label,'00')||strcmp(s.label,'99')
        sind=[];
    else
        sind = find(i==s.index);
    end
    if i(1)>3, jj = 1;else jj = i(1)-1;end
    ps = size(i(1)-jj:i(end),2);
    if (size(MC.D2,2)>0) %case 2D      
        for siz=1:size(MC.D2,2)
            figure(MC.D2(siz));
            hold on
            dat = get(MC.D2(siz),'UserData');
            ux = get(dat.traj,'XData');
            uy = get(dat.traj,'YData');
            s1=eval(dat.label{1});
            s2=eval(dat.label{2});
            if ~isfield(gds,'combineddata') || ~iscell(gds.combineddata)
                gds.combineddata = {};
                gds.combineddata1 = {};
            end
            if size(s2,1) == 1
                dook = 1;
            else
                dook = 0;
                if ~isfield(gds,'combineddata')
                    gds.combineddata{siz} = [];
                end
                if length(gds.combineddata) < siz
                    gds.combineddata{siz} = [];
                end
                uy = gds.combineddata{siz};
                gds.combineddata{siz} = [uy s2];
            end
            if size(s1,1) == 1
                dook1 = 1;
            else
                dook1 = 0;
                if ~isfield(gds,'combineddata1')
                    gds.combineddata1{1} = [];
                end
                if i(1) == 1
                    gds.combineddata1{1} = [];
                end
                ux = gds.combineddata1{siz};
                gds.combineddata1{siz} = [ux s1];
            end
            figure(MC.D2(siz));
            hold on
            if ~dook1 && ~dook
                for ii = 1:size(s1,1)
                    t1 = gds.combineddata1{siz};
                    t2 = gds.combineddata{siz};
                    line(t1(ii,:),t2(ii,:));
                    drawnow;
                end
            elseif ~dook1                
                for ii = 1:size(s1,1)
                    t1 = gds.combineddata1{siz};
                    line(t1(ii,end-length([uy s2])+1:end),[uy s2]);
                    gds.combineddata{siz} = [uy s2];
                    drawnow;
                end
            elseif ~dook     
                for ii = 1:size(s2,1)
                    t2 = gds.combineddata{siz};
                    line([ux s1],t2(ii,end-length([ux s1])+1:end));
                    gds.combineddata1{siz} = [ux s1];
                    drawnow
                end
            else
                line([ux s1],[uy s2]);
                    gds.combineddata1{siz} = [ux s1];
                    gds.combineddata{siz} = [uy s2];
                drawnow;
            end
            if ~isempty(sind) %case singularity detected  
                uxs = get(dat.trajs,'XData');
                uys = get(dat.trajs,'YData');
                s3 = [];
                s4 = [];
                tmp = dat.labels{1};
                for ij=1:size(tmp,1)
                    s3=[s3;eval(tmp(ij,:))];
                end
                tmp = dat.labels{2};
                for ij=1:size(tmp,1)
                    s4=[s4;eval(tmp(ij,:))];
                end
                line([uxs s3],[uys s4],'Linestyle', 'none' , 'Marker' , '*','Color','r');
                d = axis;
                skew = 0.01*[d(2)-d(1) d(4)-d(3)];
                if length(s3) < length(s4)
                    s3 = s3 * ones(size(s4));
                elseif length(s3) > length(s4)
                    s4 = s4 * ones(size(s3));
                end                    
                text(s3+skew(1), s4+skew(2), s.label );
            end    
            drawnow;   
            hold off
        end   
    end
    if i(1) == 1
        gds.combined3data1 = {};
        gds.combined3data2 = {};
        gds.combined3data3 = {};
    end
    if size(MC.D3,2)>0%case 3D
        for siz = 1:size(MC.D3,2)               
            figure(MC.D3(siz));
            hold on
            dat = get(MC.D3(siz),'UserData');
            ux = get(dat.traj,'XData');
            uy = get(dat.traj,'YData');
            uz = get(dat.traj,'ZData');
                s1 = eval(dat.label{1});
                s2 = eval(dat.label{2});
                s3 = eval(dat.label{3});
            if ~isfield(gds,'combined3data1') || ~iscell(gds.combined3data1)
                gds.combined3data3 = {};
                gds.combined3data2 = {};
                gds.combined3data1 = {};
            end
            if size(s1,1) == 1
                dook1 = 1;
            else
                dook1 = 0;
                if ~isfield(gds,'combined3data1')
                    gds.combined3data1{1} = [];
                end
                if i(1) == 1
                    gds.combined3data1{1} = [];
                end
                ux = gds.combined3data1{siz};
                gds.combined3data1{siz} = [ux s1];
            end
            if size(s2,1) == 1
                dook2 = 1;
            else
                dook2 = 0;
                if ~isfield(gds,'combined3data2')
                    gds.combined3data2{siz} = [];
                end
                if length(gds.combined3data2) < siz
                    gds.combined3data2{siz} = [];
                end
                uy = gds.combined3data2{siz};
                gds.combined3data2{siz} = [uy s2];
            end
            if size(s3,1) == 1
                dook3 = 1;
            else
                dook3 = 0;
                if ~isfield(gds,'combined3data3')
                    gds.combined3data3{siz} = [];
                end
                if length(gds.combined3data3) < siz
                    gds.combined3data3{siz} = [];
                end
                uz = gds.combined3data3{siz};
                gds.combined3data3{siz} = [uz s3];
            end
            figure(MC.D3(siz));
            hold on
            if ~dook1 && ~dook2 && ~dook3
                for ii = 1:size(s1,1)
                    t1 = gds.combined3data1{siz};
                    t2 = gds.combined3data2{siz};
                    t3 = gds.combined3data3{siz};
                    line(t1(ii,:),t2(ii,:),t3(ii,:));
                    drawnow;
                end
            elseif ~dook1 
                if ~dook2
                    for ii = 1:size(s1,1)
                        t1 = gds.combined3data1{siz};
                        t2 = gds.combined3data2{siz};
                        line(t1(ii,end-length([uz s3])+1:end),t2(ii,end-length([uz s3])+1:end),[uz s3]);
                        gds.combined3data3{siz} = [uz s3];
                        drawnow;
                    end
                elseif ~dook3
                    for ii = 1:size(s1,1)
                        t1 = gds.combined3data1{siz};
                        t3 = gds.combined3data3{siz};
                        line(t1(ii,end-length([uy s2])+1:end),[uy s2],t3(ii,end-length([uy s2])+1:end));
                        gds.combined3data2{siz} = [uy s2];
                        drawnow;
                    end
                else
                    for ii = 1:size(s1,1)
                        t1 = gds.combined3data1{siz};
                        line(t1(ii,end-length([uy s2])+1:end),[uy s2],[uz s3]);
                        gds.combined3data2{siz} = [uy s2];
                        gds.combined3data3{siz} = [uz s3];
                        drawnow;
                    end
                end
            elseif ~dook2   
                if ~dook3
                    for ii = 1:size(s2,1)
                        t2 = gds.combined3data2{siz};
                        t3 = gds.combined3data3{siz};
                        line([ux s1],t2(ii,end-length([ux s1])+1:end),t3(ii,end-length([ux s1])+1:end));
                        gds.combined3data1{siz} = [ux s1];
                        drawnow
                    end
                else
                    for ii = 1:size(s2,1)
                        t2 = gds.combined3data2{siz};
                        line([ux s1],t2(ii,end-length([ux s1])+1:end),[uz s3]);
                        gds.combined3data1{siz} = [ux s1];
                        gds.combined3data3{siz} = [uz s3];
                        drawnow
                    end
                end
            else
                if ~dook3
                    for ii = 1:size(s3,1)
                        t3 = gds.combined3data3{siz};
                        line([ux s1],[uy s2],t3(ii,end-length([ux s1])+1:end));
                        gds.combined3data1{siz} = [ux s1];
                        gds.combined3data2{siz} = [uy s2];
                        drawnow;
                    end
                else
                    line([ux s1],[uy s2],[uz s3]);
                    gds.combined3data1{siz} = [ux s1];
                    gds.combined3data2{siz} = [uy s2];
                    gds.combined3data3{siz} = [uz s3];
                    drawnow;
                end
            end
            if ~isempty(sind) %case singularity detected  
                uxs = get(dat.trajs,'XData');
                uys = get(dat.trajs,'YData');
                uzs = get(dat.trajs,'ZData');
                s4 = []; s5 = []; s6 = [];
                s1 = eval(dat.label{1});s2  = eval(dat.label{2});  s3 = eval(dat.label{3});
                tmp = dat.labels{1};
                for ij=1:size(tmp,1)
                    s4=[s4;eval(tmp(ij,:))];
                end
                tmp = dat.labels{2};
                for ij=1:size(tmp,1)
                    s5=[s5;eval(tmp(ij,:))];
                end
                tmp = dat.labels{3};
                for ij=1:size(tmp,1)
                    s6=[s6;eval(tmp(ij,:))];
                end
                mx = max([length(s4),length(s5),length(s6)]);
                if mx > length(s4)
                    s4 = s4 * ones(mx,1);
                end             
                if mx > length(s5)
                    s5 = s5 * ones(mx,1);
                end       
                if length(s6) < mx
                    s6 = s6 * ones(mx,1);
                end                 
                line([uxs s4],[uys s5],[uzs s6],'Linestyle', 'none' , 'Marker' , '*','Color','r');
                figure(MC.D3(siz));
                d = axis;
                skew = 0.01*[d(2)-d(1) d(4)-d(3) d(6)-d(5)]; 
                text(s4+skew(1),s5+skew(2), s6+skew(3),s.label );
            end    
            drawnow;   
            hold off
        end   
    end      
    %numeric window is open
    if ~isempty(MC.numeric_fig)
        for d=1:ndim
            parameters(d) = gds.parameters{d,2};
        end   
        j = gds.options.ActiveParams(1,1);
        parameters(j) = xout((gds.dim)+1,i(end));
        dat = get(MC.numeric_fig,'Userdata');
        for k=1:size(dat.label,2)
                eval(dat.label{k});  
        end
    end
%------------------------------------------------------    
function load_point(index,x,string,file,varargin)
global gds eds
    load(file);
    for i=1:gds.dim
        gds.coordinates{i,2} = x(i,index);
    end
    dim = size(eds.P0,1);
    for i=1:dim
        gds.parameters{i,2} = eds.P0(i,1);
    end    
    j = eds.ActiveParams(1,1);
    gds.parameters{j,2} = x((gds.dim)+1,index);
    if strcmp(string,'save')
        num = varargin{1};
        save(file,'x','v','s','h','f','num','eds','cds','ctype','point');
    end
 %-----------------------------------------------------
 function singularities
 global gds 
    if size(gds.options.IgnoreSingularity,2)==3
        gds.options = contset(gds.options,'Singularities',0);
    else
        gds.options = contset(gds.options,'Singularities',1);
    end     
%------------------------------------------------------    
function start_cont(varargin)
global gds path_sys cds eds MC
    if ~isempty(varargin)
        file = fullfile(path_sys,gds.system,gds.diagram,strcat(gds.curve.new,'.mat'));
        if exist(file,'file')
            load(file);
        else
            errordlg('It is not possible to extend the current curve!','Error extend EP');
            return 
        end  
        [x,v,s,h,f] = cont(x,v,s,h,f,cds);
    else
        if (size(gds.options.ActiveParams,2)~=1)
            errordlg('You have to select exactly one free parameter');
            return
        end
        x0  = vertcat(gds.coordinates{:,2});
        par = vertcat(gds.parameters{:,2});
        systemhandle = str2func(gds.system);
        switch gds.point
        case 'BP '
            if ~isempty(gds.curve.old)
                file = strcat(gds.curve.old,'.mat');
                file = fullfile(path_sys,gds.system,gds.diagram,file);
                if exist(file,'file')
                    load(file);
                else
                    errordlg('It is not possible to start from this point. There isn''t enough information available!','Error start cont EP');
                    return 
                end  
            else
                errordlg('It is not possible to start from this point. There isn''not enough information available!','Error start cont EP');
                return
            end                   
            [x0,v0] = init_BP_EP(systemhandle,x0,par,s(num),gds.amplitude);
        case 'LP '
            [x0,v0] = init_LP_EP(systemhandle,x0,par,gds.options.ActiveParams);
        case 'H '
            [x0,v0] = init_H_EP(systemhandle,x0,par,gds.options.ActiveParams);
        otherwise
            if isempty(gds.curve.old)
                gds.curve.old = gds.curve.new;
            end
            [x0,v0] = init_EP_EP(systemhandle,x0,par,gds.options.ActiveParams);
        end
        [x,v,s,h,f] = cont(@equilibrium,x0,v0,gds.options);
    end
    status = mkdir(path_sys,gds.system);
    dir = fullfile(path_sys,gds.system);
    status = mkdir(dir,gds.diagram);
    file = fullfile(path_sys,gds.system,gds.diagram,gds.curve.new);
    point = gds.point; ctype = gds.type; 
    if ~exist('num','var'),num = 1;end
    Adapt = contget(gds.options,'Adapt',3);
    cds.options = contset(cds.options,'Adapt',Adapt);
    if ~isempty(x)
        save(file,'x','v','s','h','f','cds','eds','ctype','point','num');                      
    end


%---------------------------------------------------------------------------
function j = start_monitor(handles,j)
global gds path_sys;
color = [1 1 1];
j = j-2;
pos  = [10 j 35 1.80];user.pos=pos;user.num=0;
stat = uicontrol(handles.figuur,'Style','text','String','Monitor Singularities','Tag','singularities','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',pos,'UserData',user);
%BP=1;H=1;LP=3
for k = 1:3
    tag1 = strcat('text',num2str(k));           
    switch k
    case 1
        string = 'Branching';
    case 2
        if (gds.dim < 2),
            continue;
        end
        string = 'Hopf';
    case 3
        string = 'Limit point';       
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
function j = start_testfunctions(handles,j)
global gds ;
color = [1 1 1];
j = j-2;
pos = [5 j 38 1.80]; user.num = 0; user.pos = pos;
string = strcat('Testfunctions (',gds.type,')');
stat = uicontrol(handles.numeric_fig,'Style','text','String',string,'Tag','testfunctions','BackGroundColor',color,'Units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',pos,'UserData',user);
%BP=1;H=2;LP=3
for k = 1:3
    j = j-2;
    post = [2 j 18 1.80]; user.pos = post;
    pose = [20 j 24 1.80];     

    switch k
    case 1
        tag = 'BP';
    case 2
        if (gds.dim< 2),
            continue;
        end
        tag = 'H';
    case 3
        tag = 'LP';
    otherwise
        errordlg('Is not possible');
    end
    stat1 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String',tag,'BackGroundColor',color,'Units','characters','fontname','FixedWidth','fontsize',12);
    stat2 = uicontrol(handles.numeric_fig,'Style','text','HorizontalAlignment','left','String','','Tag',tag,'BackGroundColor',color,'Units','characters','fontname','FixedWidth','fontsize',12);
    set(stat1,'Position',post,'UserData',user); user.pos = pose;
    set(stat2,'Position',pose,'UserData',user);
end
guidata(handles.numeric_fig,handles);

%---------------------------------------------------------------------------------------
function j = start_amplitude(handles,j)
%enter field for amplitude
global gds path_sys;
color = [1 1 1];
j = j-2;
pos  = [10 j 35 1.80]; user.num = 0; user.pos = pos;
stat = uicontrol(handles.figuur,'Style','text','String','Switch Data','Tag','discretization','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',pos,'UserData',user);
j = j-2;
post = [2 j 18 1.80];user.pos=post;
pose = [20 j 25 1.80];
if ~isfield('gds','amplitude')
    gds.amplitude = 1e-6;
    file = fullfile(path_sys,gds.system);
    save(file,'gds');
end
stat = uicontrol(handles.figuur,'Style','text','HorizontalAlignment','left','String','amplitude','Tag','text_amplitude','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
set(stat,'Position',post,'UserData',user);
edit = uicontrol(handles.figuur,'Style','edit','HorizontalAlignment','left','String',gds.amplitude,'Tag','amplitude','BackGroundColor',color,'units','characters','fontname','FixedWidth','fontsize',12);
user.pos = pose; set(edit,'Position',pose,'Callback','set_option','UserData',user);
guidata(handles.figuur,handles);
