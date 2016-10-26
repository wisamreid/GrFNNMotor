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

  
global gds sys ud MC driver_window calculation_progress string_O npoints sp ncross integDOstr1 integDOstr2
status = 0;
nt = length(t);
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
        
        
        % Poincare map analysis
        k=11;     
        if (gds.poincare_do>0)
            % Poincare map by secant to both direction                  
            X = y(:,1)';
            XIK = X;
            t = t(1);
            tik = t;
            rr = feval(gds.poincare_eq,t,X,gds.parameters{:,2});
            rr1=rr;XK=X;TK=t;
            XL=ud.y(newi-1,:);
            TL=ud.t(newi-1);
            if sign(rr)~=sign(ud.poincare_cur)  
            % If the testfunction has changed sign
                if (gds.integrator.nr_crossection == 0) || (sign(rr) == sign(gds.integrator.nr_crossection))
                % Then this intersection is supposed to be shown, because
                % all intersections are asked or this particular direction
                    ncross=0;
                    k = 0;   
                    while (abs(rr)>1e-3) && (k<=10)
                        k = k+1;
                        X = feval(gds.integrator_fun,TK,XK,gds.parameters{:,2});     
                        t = TK;
                        df_plane=feval(gds.poincare_eq,t,X, gds.parameters{:,2})-ud.poincare_dif;
                        TK1 = TK - rr/df_plane;
                        Poinc_options = odeset('RelTol',gds.integrator.options.RelTol,'AbsTol',gds.integrator.options.AbsTol,'MaxStep',gds.integrator.options.MaxStep,...
                            'Refine',2,'InitialStep',TK1 - TL);  
                        [TT,XT]=feval(gds.integrator.method,gds.integrator_fun,[TL TK1],XL,Poinc_options,gds.parameters{:,2});
                        kXT = size(XT);
                        t = TT(kXT(1));
                        X = XT(kXT(1),:);
                        rr = feval(gds.poincare_eq,t,X, gds.parameters{:,2});
                        XK = X;
                        TK = t;
                        if (rr*ud.poincare_cur)>0
                            XL = X;
                            TL = t;
                        end
                    end 
                end
            end
            if k<11 
                %insert value of poincaré
                 y(:,1)=X;
                calculation_progress=2;
                ud.poincare_cur = rr1; 
                ud.poincare_nmbr = ud.poincare_nmbr + 1;
                sp(end+1,1).index = 1;
                if length(sp) > 1
                    sp(end,1).index = sp(end-1,1).index+1;
                end
                sp(end,1).data.t = TK;
                sp(end,1).data.x = XK;
                sp(end,1).msg ='Poincaré';
                sp(end,1).label ='P ';
                nwind = size(MC.D2,2); 
                    for i=1:nwind 
                        ff = get(MC.D2(i),'UserData');
                        Xp = eval(ff.labels{1});
                        Yp = eval(ff.labels{2});
                        Xp = XK(1);
                        Yp = XK(2);
                        if ud.poincare_nmbr==1 
                            set(ff.poinc_trj,'Xdata',Xp,'Ydata',Yp,'color','r');
                        else 
                            ux=get(ff.poinc_trj,'XData');
                            uy=get(ff.poinc_trj,'YData');
                            set(ff.poinc_trj,'Xdata',[ux Xp'],'Ydata',[uy Yp'],'color','r');
                        end 
                        set(ff.poinc_trj,'Visible','on');
                        drawnow; 
                    end
                    nwind = size(MC.D3,2); 
                    for i=1:nwind 
                        ff = get(MC.D3(i),'UserData');
                        Xp = eval(ff.labels{1});
                        Yp = eval(ff.labels{2});
                        Zp = eval(ff.labels{3});
                        Xp = XK(1);
                        Yp = XK(2);
                        Zp = XK(3);
                        if ud.poincare_nmbr==1 
                            set(ff.poinc_trj,'Xdata',Xp,'Ydata',Yp,'Zdata',Zp,'color','r');
                        else 
                            ux=get(ff.poinc_trj,'XData');
                            uy=get(ff.poinc_trj,'YData');
                            uz=get(ff.poinc_trj,'ZData');
                            set(ff.poinc_trj,'Xdata',[ux Xp'],'Ydata',[uy Yp'],'Zdata',[uz Zp'],'color','r');
                        end 
                        set(ff.poinc_trj,'Visible','on');
                        drawnow; 
                    end
                    if ~isempty(MC.numeric_fig)%numeric window is open
                        ff = get(MC.numeric_fig,'Userdata');
                        for k=1:size(dat.label,2)
                            eval(ff.labels{k});  
                        end
                    end 
%                     PRS(2,3,MC.mainwindow.mstatus);
            end
            X = XIK;
            t = tik;
            ad.poincare_nmbr=ud.poincare_nmbr;
            if gds.poincare_do>0
                if ~isnan(X(1))
                    ud.poincare_cur=feval(gds.poincare_eq,t,X,gds.parameters{:,2});
                end
            end 
            ad.poincare_cur=ud.poincare_cur;
        end
        PRS(sp(end).index,newi,MC.mainwindow.mstatus);
        ud = ad; 
        
    else % no output at this time
        return;
    end  
       
else
    switch(flag)
    case 'init'                                % odephas2(tspan,y0,'init')
        nt=1;    
        npoints=sys.gui.plot_points;
        ncross=0;
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
        if ~isfield(gds,'poincare_do'),gds.poincare_do=0;end
        if gds.poincare_do>0
            Y = y(:,1)';
            t = t(1);
            ud.poincare_cur=feval(gds.poincare_eq,t,Y, gds.parameters{:,2});
            Y(:) = 0;
            ud.poincare_dif=feval(gds.poincare_eq,t,Y, gds.parameters{:,2});
            ud.poincare_nmbr = 0;
        else
%             errordlg('Please enter a function!');
%             return;
            ud.poincare_cur=[];ud.poincare_dif=[];ud.poincare_nmbr=[];
        end
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
                ud.labels{1}='empt';ud.labels{2}='empt';                
                for k=1:2
                    if ud.p2(1,k)<0        
                        pl2=-(ud.p2(1,k));
                        ud.label{k} = sprintf('gds.parameters{%d,2}*ones(nt,1)',pl2);
                        ud.labels{k} = sprintf('gds.parameters{%d,2}',pl2);
                    end                 
                    if ud.p2(1,k)==0
                        ud.label{k} ='ud.t(1:newi,1)';
                        ud.labels{k} ='tout(siz)';
                    end
                    if strcmp(ud.label{k},'empt')==1       
                        ud.label{k} = sprintf('ud.y(1:newi,%d)',ud.p2(1,k));
                    end
                    if strcmp(ud.labels{k},'empt')==1       
%                         ud.labels{k} = sprintf('yout(%d,siz)',ud.p2(1,k));
                         ud.labels{k} =  sprintf('ud.y(1:newi,%d)',ud.p2(1,k));
                    end
                end
                Xp = eval(ud.label{1});
                Yp = eval(ud.label{2});
                ud.traj = line(Xp,Yp,'Visible','on');
                ud.poinc_trj = line(Xp,Yp,'Marker','.','LineStyle','none','MarkerSize',5,'Color',[1 1 1]);
                set(ud.poinc_trj,'Xdata',[],'Ydata',[]);
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
                ud.labels{1}='empt';ud.labels{2}='empt';ud.labels{3}= 'empt';
                for k=1:3
                    if ud.p3(1,k)<0        
                        pl3=-(ud.p3(1,k));
                        ud.label{k}=sprintf('gds.parameters{%d,2}*ones(nt,1)',pl3);
                    end   
                    if ud.p3(1,k)==0
                        ud.label{k}='ud.t(1:newi,1)';
                        ud.labels{k} ='t(1)';
                    end
                    if strcmp(ud.label{k},'empt')==1       
                        ud.label{k}=sprintf('ud.y(1:newi,%d)',ud.p3(1,k));                
                    end
                    if strcmp(ud.labels{k},'empt')==1       
                        ud.labels{k} = sprintf('X(%d)',ud.p3(1,k));
                    end
                end   
                Xp = eval(ud.label{1});
                Yp = eval(ud.label{2});
                Zp = eval(ud.label{3});
                ud.traj = line(Xp,Yp,Zp);
                set(ud.traj,'XData',[],'YData',[],'ZData',[]);
                ud.poinc_trj = line(Xp,Yp,'Marker','.','LineStyle','none','MarkerSize',5,'Color',[1 1 1]);
                set(ud.poinc_trj,'Xdata',[],'Ydata',[],'ZData',[]);
                xlabel(gds.plot3(nr,1).label); ylabel(gds.plot3(nr,2).label);zlabel(gds.plot3(nr,3).label);
                set(MC.D3(siz),'UserData',ud); 
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
            if ~(size(MC.D3,2)==1),d=vertcat(d{:,1});
            end
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
                for k=1:size(y,2)
                    for d=1:gds.dim
                        ud.label{num}=sprintf('set(MC.numeric_handles.coord%d,''String'',num2str(y(%d,%d),''%%.8g''))',d,d,k);
                        ud.labels{num}=sprintf('set(MC.numeric_handles.coord%d,''String'',num2str(y(%d,end),''%%.8g''))',d,d);
                        num = num+1;
                    end
                end
            end
            if isfield(MC.numeric_handles,'time')
                ud.label{num}=sprintf('set(MC.numeric_handles.time,''String'',num2str(t(%d,1),''%%.10g''))',k); 
                ud.labels{num}=sprintf('set(MC.numeric_handles.time,''String'',num2str(t(end),''%%.10g''))');
                num = num+1;
            end
            if isfield(MC.numeric_handles,'param1')
                for d=1:size(gds.parameters,1)
                   ud.label{num}=sprintf('set(MC.numeric_handles.param%d,''String'',gds.parameters{%d,2})',d,d);
                   ud.labels{num}=ud.label{num};
                   num = num+1;
                end 
            end
            set(MC.numeric_fig,'Userdata',ud);
        end
        switch gds.integrator.method
        case 'ode45'
            integDOstr1='h=  evalin(''caller'',''h'');f=  evalin(''caller'',''f'');';
            integDOstr2='[tout,yout,iout,vnew,found]=integDOzero(@mntrp45,v1,t1,y1,tend,yend,gds.integrator.tspan(1),h,f);';
        case 'ode23'
            integDOstr1='h=  evalin(''caller'',''h'');f=  evalin(''caller'',''f'');';
            integDOstr2='[tout,yout,iout,vnew,found]=integDOzero(@mntrp23,v1,t1,y1,tend,yend,gds.integrator.tspan(1),h,f);';
        case 'ode113'
            integDOstr1='klast=  evalin(''caller'',''klast'');phi=  evalin(''caller'',''phi'');psi=  evalin(''caller'',''psi'');';
            integDOstr2='[tout,yout,iout,vnew,found]=integDOzero(@mntrp113,v1,t1,y1,tend,yend,gds.integrator.tspan(1),klast,phi,psi);';
        case 'ode15s'
            integDOstr1='h=  evalin(''caller'',''h'');dif=  evalin(''caller'',''dif'');k=  evalin(''caller'',''k'');';
            integDOstr2='[tout,yout,iout,vnew,found]=integDOzero(@mntrp15s,v1,t1,y1,tend,yend,gds.integrator.tspan(1),h,dif,k);';
        case 'ode23s'
            integDOstr1='h=  evalin(''caller'',''h'');k1=  evalin(''caller'',''k1'');k2=  evalin(''caller'',''k2'');';
            integDOstr2='[tout,yout,iout,vnew,found]=integDOzero(@mntrp23s,v1,t1,y1,tend,yend,gds.integrator.tspan(1),h,k1,k2);';
        case 'ode23t'
            integDOstr1='h=  evalin(''caller'',''h'');z=  evalin(''caller'',''z'');znew=  evalin(''caller'',''znew'');';
            integDOstr2='[tout,yout,iout,vnew,found]=integDOzero(@mntrp23t,v1,t1,y1,tend,yend,gds.integrator.tspan(1),h,z,znew);';
        case 'ode23tb'
            integDOstr1='t2=  evalin(''caller'',''t2'');y2=  evalin(''caller'',''y2'');';
            integDOstr2='[tout,yout,iout,vnew,found]=integDOzero(@mntrp23tb,v1,t1,y1,tend,yend,gds.integrator.tspan(1),t2,y2);';
%         case 'ode78'
%             str = strcat('ode78',str);
%         case 'ode87'
%             str = strcat('ode87',str);
        end  
        set(MC.mainwindow.mstatus,'String','computing');
        % The STOP button.        
        stop;
        
        
             
        

    case 'done'         % integ_plot([],[],'done')   
        newi = ud.i; 
        
        if size(MC.D2,2)>0
            for siz=1:size(MC.D2,2)
                dat = get(MC.D2(siz),'UserData');
                ux=get(dat.traj,'XData');
                uy=get(dat.traj,'YData');
                s1=eval(dat.label{1});
                nt=length(s1);
                s2=eval(dat.label{2});
                set(dat.traj,'XData',[ux s1'],'YData',[uy s2']);
                drawnow;   
            end
        end
        if size(MC.D3,2)>0
            for siz=1:size(MC.D3,1)
                dat = get(MC.D3(siz),'UserData');
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
            ud = rmfield(ud,{'t','y','i','poincare_cur','poincare_dif','poincare_nmbr','poinc_trj'});
            for k=1:size(MC.D2,2)
                ud(k).p2=[];ud(k).label=[];ud(k).traj=[];ud(k).labels=[];ud(k).trajs=[];
                set(MC.D2(k),'UserData',ud(k));
            end
        elseif size(MC.D3,2)>0
            ud = get(MC.D3,'UserData');
            if (size(MC.D3,2)~=1),ud=vertcat(ud{:,1});end
            ud = rmfield(ud,{'t','y','i','poincare_cur','poincare_dif','poincare_nmbr','poinc_trj'});
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
