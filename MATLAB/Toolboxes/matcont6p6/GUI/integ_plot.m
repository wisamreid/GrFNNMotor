function status = integ_plot(t,y,flag,varargin)
%    integ_plot  2and3-D phase plane ODE output function.
%   When the function is passed to an ODE solver as the 'OutputFcn'
%   property, i.e. options = odeset('OutputFcn',@odephas2), the solver
%   calls ODEPHAS2(T,Y,'') after every timestep.     
%   At the start of integration, a solver calls INTEG_PLOT(TSPAN,Y0,'init') to
%   initialize the output function.  After each integration step to new time
%   point T with solution vector Y the solver calls STATUS = INTEG_PLOT(T,Y,'').
%   If the solver's 'Refine' property is greater than one (see ODESET), then
%   T is a column vector containing all new output times and Y is an array
%   comprised of corresponding column vectors.  The STATUS return value is 1
%   if the STOP button has been pressed and 0 otherwise.  When the
%   integration is complete, the solver calls INTEG_PLOT([],[],'done').

  
global gds sys ud MC driver_window calculation_progress string_O npoints;
status = 0;
if isempty(flag) % odephas2(t,y) [v5 syntax] or odephas2(t,y,'')
    nt = length(t);
    chunk = max(128,nt);
    [rows,cols] = size(ud.y);
    oldi = ud.i;   
    newi = oldi + nt;
    if newi > rows
         ud.t = [ud.t; zeros(chunk,1)];
         ud.y = [ud.y; zeros(chunk,cols)];
    end     
    ud.t(oldi+1:newi) = t;
    ud.y(oldi+1:newi,:) = y.';
    ud.i = newi;
    if calculation_progress==0
       status = 1;
       return;
    else 
       status = 0;
    end;
    if calculation_progress==2
        set(MC.mainwindow.mstatus,'String','pause');
        waitfor(driver_window,'UserData',0);
        set(MC.mainwindow.mstatus,'String','computing');
    end;
    if (fix(ud.i/npoints)>0)
        ad = ud;
        ad.t = zeros(chunk,1);
        ad.y = zeros(chunk,gds.dim);
        ad.t(1) = t(end);
        ad.y(1,:) = y(:,end)';
        ad.i = 1;
        if size(MC.D2,2)>0
            for siz=1:size(MC.D2,2)
                dat = get(MC.D2(siz),'UserData');
                ux = get(dat.traj,'XData');
                uy = get(dat.traj,'YData');
                s1=eval(dat.label{1});
                s2=eval(dat.label{2});
                set(dat.traj,'XData',[ux s1'],'YData',[uy s2']);
                set(dat.traj,'Visible','on');
                drawnow;   
            end
        end
        if size(MC.D3,2)>0
            for siz=1:size(MC.D3,2)
                dat = get(MC.D3(siz),'UserData');
                ux = get(dat.traj,'XData');
                uy = get(dat.traj,'YData');
                uz = get(dat.traj,'ZData');
                s1 = eval(dat.label{1});
                s2 = eval(dat.label{2});
                s3 = eval(dat.label{3});
                set(dat.traj,'XData',[ux s1'],'YData',[uy s2'],'ZData', [uz s3']);
                drawnow;   
            end
        end
        if ~isempty(MC.numeric_fig)%numeric window is open
            y=ud.y(newi,:)';
            t=ud.t(newi);          
            dat = get(MC.numeric_fig,'Userdata');
            for k=1:size(dat.label,2)
                eval(dat.label{k});  
            end
        end   
        ud = ad;                  
    else % no output at this time
        return;
    end  
       
else
    switch(flag)
    case 'init'                           % odephas2(tspan,y0,'init')
        npoints=sys.gui.plot_points;
        string_O = feval(gds.gui_load_draw);
        ud = [];
        ud.nr = 1;
        cols = length(y);
        ud.t = zeros(128,1);
        ud.y = zeros(128,gds.dim);
        ud.i = 1;
        newi = 1;
        oldi = 1;
        ud.t(1) = t(1);
        ud.y(1,:) = y(1:gds.dim).';
        % Rather than redraw all data at every timestep, we will simply move
        % the last line segment along, not erasing it.
        
        % The STOP button.        
       for i=1:gds.dim,x(i,1)=0;end;    
        if size(MC.D2,2)>0
            ud.p2 = [];
            for siz=1:size(MC.D2,2)
                figure(MC.D2(siz));
                dat = get(MC.D2(siz),'UserData');
                ud.nr = dat.nr;nr=dat.nr;
                if isempty(strmatch(char(gds.plot2(nr,1).type),string_O,'exact'))|isempty(strmatch(char(gds.plot2(nr,2).type),string_O,'exact'))
                    ud.p2=[inf inf];
                else
                    ud.p2 =feval(gds.gui_make_ready,gds.plot2(nr,:),x);           
                end
                ud.label{1}='empt';ud.label{2}='empt';
                for k=1:2
                    if ud.p2(1,k)<0        
                        pl2=-(ud.p2(1,k));
                        ud.label{k} = sprintf('gds.parameters{%d,2}*ones(1,1)',pl2);                      
                    end                 
                    if ud.p2(1,k)==0
                        ud.label{k} ='ud.t(1:newi,1)';
                    end
                    if strcmp(ud.label{k},'empt')==1   
                        ud.label{k} = sprintf('ud.y(1:newi,%d)',ud.p2(1,k));
                    end
                end
                Xp = eval(ud.label{1});
                Yp = eval(ud.label{2});
                ud.traj = line(Xp,Yp);
                set(ud.traj,'XData',[],'YData',[]);
                xlabel(gds.plot2(nr,1).label);ylabel(gds.plot2(nr,2).label);
                set(MC.D2(siz),'UserData',ud); 
            end
        end
        if size(MC.D3,2)>0
            ud.p3 = [];
            for siz=1:size(MC.D3,2)
                figure(MC.D3(siz));
                dat = get(MC.D3(siz),'UserData');
                ud.nr = dat.nr;nr=dat.nr;
                if isempty(strmatch(char(gds.plot3(nr,1).type),string_O,'exact'))|isempty(strmatch(char(gds.plot3(nr,2).type),string_O,'exact'))|isempty(strmatch(char(gds.plot3(nr,3).type),string_O,'exact'))
                    ud.p3=[inf inf inf];
                else
                    ud.p3 =feval(gds.gui_make_ready,gds.plot3(nr,:),x);           
                end
                ud.label{1}='empt';ud.label{2}='empt';ud.label{3} = 'empt';
                for k=1:3
                    if ud.p3(1,k)<0        
                        pl3=-(ud.p3(1,k));
                        ud.label{k}=sprintf('gds.parameters{%d,2}*ones(nt,1)',pl3);
                    end   
                    if ud.p3(1,k)==0
                        ud.label{k}='ud.t(1:newi,1)';
                    end
                    if strcmp(ud.label{k},'empt')==1       
                        ud.label{k}=sprintf('ud.y(1:newi,%d)',ud.p3(1,k));                
                    end
                end   
                Xp = eval(ud.label{1});
                Yp = eval(ud.label{2});
                Zp = eval(ud.label{3});
                ud.traj = line(Xp,Yp,Zp);
                set(ud.traj,'XData',[],'YData',[],'ZData',[]);
                set(MC.D3(siz),'UserData',ud);
                xlabel(gds.plot3(nr,1).label); ylabel(gds.plot3(nr,2).label);zlabel(gds.plot3(nr,3).label);
            end
        end
        if size(MC.D2,2)>0 %2D window open
            d=get(MC.D2,'UserData');
            if ~(size(MC.D2,2)==1),d=vertcat(d{:,1});end
            da=vertcat(d.p2);
            hel=find(da(:,1)==inf);
            MC.D2(hel)=[];
            if isempty(MC.D2),MC.D2=[];end
        end
        if size(MC.D3,2)>0%3D window open
            d=get(MC.D3,'UserData');   
            if ~(size(MC.D3,2)==1),d=vertcat(d{:,1});end
            da=vertcat(d.p3);        
            hel=find(da(:,1)==inf);MC.D3(hel)=[];
            if isempty(MC.D3),MC.D3=[];end
        end
        if ~isempty(MC.numeric_fig)
            ud.label{1}='';
            num = 1;
            dat = get(MC.numeric_fig,'UserData');
            ud.pos = dat.pos;
            if isfield(MC.numeric_handles,'coord1')
                    for d=1:gds.dim
                        ud.label{num}=sprintf('set(MC.numeric_handles.coord%d,''String'',num2str(y(%d,end),''%%.10g''))',d,d);
                        num = num+1;
                    end
            end
            if isfield(MC.numeric_handles,'time')
                ud.label{num}='set(MC.numeric_handles.time,''String'',num2str(t(end,1)))';          
                num = num+1;
            end
            if isfield(MC.numeric_handles,'param1')
                for d=1:size(gds.parameters,1)
                   ud.label{num}=sprintf('set(MC.numeric_handles.param%d,''String'',gds.parameters{%d,2})',d,d);
                   num = num+1;
                end 
            end            
            set(MC.numeric_fig,'Userdata',ud);
        end
        set(MC.mainwindow.mstatus,'String','computing');
        stop;

    case 'done'         % integ_plot([],[],'done')  
        newi = ud.i;
        if size(MC.D2,2)>0
            for siz=1:size(MC.D2,2)
                dat = get(MC.D2(siz),'UserData');
                ux=get(dat.traj,'XData');
                uy=get(dat.traj,'YData');
                s1=eval(dat.label{1});
                s2=eval(dat.label{2});
                set(dat.traj,'XData',[ux s1'],'YData',[uy s2']);
                drawnow;   
            end
        end
        if size(MC.D3,2)>0
            for siz=1:size(MC.D3,1)
                dat = get(MC.D3(siz),'UserData');
                for k=1:3
                    ux = get(dat.traj,'XData');
                    uy = get(dat.traj,'YData');
                    uz = get(dat.traj,'ZData');
                    s1 = eval(dat.label{1});
                    s2 = eval(dat.label{2});
                    s3 = eval(dat.label{3});
                    set(dat.traj,'XData',[ux s1'],'YData',[uy s2'],'ZData',[uz s3']);
                    drawnow;   
                end
            end
        end        
        if ~isempty(MC.numeric_fig)%numeric window is open
            y = ud.y(newi,:)';
            t = ud.t(newi);
            dat = get(MC.numeric_fig,'Userdata');
            for k=1:size(dat.label,2)
                eval(dat.label{k});  
            end
        end
        if size(MC.D2,2)>0
            ud = get(MC.D2,'UserData');
            if (size(MC.D2,2)~=1),ud=vertcat(ud{:,1});end
            ud = rmfield(ud,{'t','y','i'});
            for k=1:size(MC.D2,2)
                ud(k).p2=[];ud(k).label=[];ud(k).traj=[];ud(k).labels=[];ud(k).trajs=[];
                set(MC.D2(k),'UserData',ud(k));
            end
        elseif size(MC.D3,2)>0
            ud = get(MC.D3,'UserData');
            if (size(MC.D3,2)~=1),ud=vertcat(ud{:,1});end
            ud = rmfield(ud,{'t','y','i'});
            for k=1:size(MC.D3,2)
                ud(k).p3=[];ud(k).label=[];ud(k).traj=[];ud(k).labels=[];ud(k).trajs=[];
                set(MC.D3(k),'UserData',ud(k));
            end
        end                   
        set(MC.mainwindow.mstatus,'String','ready');
        if ishandle(driver_window),delete(driver_window);end
    end
end
drawnow;
