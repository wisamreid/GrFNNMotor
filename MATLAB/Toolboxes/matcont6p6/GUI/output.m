function output(numsing,xout,s,hout,fout,i)
global gds cds MC;%lds
if (i==1)
    d.nr=1;d3.nr=1;  
    d.p2=[];d3.p3=[]; 
    string=feval(gds.gui_load_draw);
    if (size(MC.D2,2)>0)%initialize case 2D
        tem=[];
        for ss=1:size(MC.D2,2) 
            figure(MC.D2(ss));
            d = get(MC.D2(ss),'UserData');
            if isempty(strmatch(char(gds.plot2(d.nr,1).type),string,'exact'))|isempty(strmatch(char(gds.plot2(d.nr,2).type),string,'exact'))
                delete(MC.D2(ss));%MC.D2(ss)=[];
                tem = [tem ss];
                continue;
            else
                d.p2 = feval(gds.gui_make_ready,gds.plot2(d.nr,:),xout(:,i));
            end
            [d.label d.labels] = feval(gds.gui_label,d.p2);
            d.traj = line(0,0);
            d.trajs = line(0,0);
            set(d.traj,'XData',[],'YData',[]);
            set(d.trajs,'XData',[],'YData',[]);
            set(MC.D2(ss),'UserData',d); 
        end
        MC.D2(tem)=[];
        if isempty(MC.D2),MC.D2=[];end
        
    end
    if (size(MC.D3,2)>0)%initialize case 3D
        tem =[];
        for ss=1:size(MC.D3,2)
            figure(MC.D3(ss));
            d3 = get(MC.D3(ss),'UserData');
            if isempty(strmatch(char(gds.plot3(d3.nr,1).type),string,'exact'))|isempty(strmatch(char(gds.plot3(d3.nr,2).type),string,'exact'))|isempty(strmatch(char(gds.plot3(d3.nr,3).type),string,'exact'))
                delete(MC.D3(ss));tem = [tem ss];
                continue;
            else
                d3.p3 = feval(gds.gui_make_ready,gds.plot3(d3.nr,:),xout(:,i));
            end
            [d3.label d3.labels] = feval(gds.gui_label,d3.p3);
            d3.traj = line(0,0,0);
            d3.trajs = line(0,0,0);
            set(d3.traj,'XData',[],'YData',[],'ZData',[]);
            set(d3.trajs,'XData',[],'YData',[],'ZData',[]);
            set(MC.D3(ss),'UserData',d3); 
        end
        MC.D3(tem)=[];
        if isempty(MC.D3),MC.D3=[];end
    end
    if ~isempty(MC.numeric_fig)
        d = get(MC.numeric_fig,'Userdata');
        d.label=feval(gds.gui_numeric_label,numsing);
        set(MC.numeric_fig,'Userdata',d);
    end
    if (size(MC.PRC,2)>0)%initialize case 2D
        for ss=size(MC.PRC,2) 
             figure(MC.PRC(ss));
            d.traj = line(0,0);
            d.trajs = line(0,0);
            set(d.traj,'XData',[],'YData',[]);
            set(d.trajs,'XData',[],'YData',[]);
            set(MC.PRC(ss),'UserData',d); 
        end
        if isempty(MC.PRC),MC.PRC=[];end
        
    end
    if (size(MC.dPRC,2)>0)%initialize case 2D
        for ss=size(MC.dPRC,2) 
             figure(MC.dPRC(ss));
            d.traj = line(0,0);
            d.trajs = line(0,0);
            set(d.traj,'XData',[],'YData',[]);
            set(d.trajs,'XData',[],'YData',[]);
            set(MC.dPRC(ss),'UserData',d); 
        end
        if isempty(MC.dPRC),MC.dPRC=[];end
        
    end
    return
end

feval(gds.gui_output,numsing,xout,s,hout,fout,i);
