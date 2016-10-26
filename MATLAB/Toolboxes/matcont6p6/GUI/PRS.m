function  npoints=PRS(index,i,status)
% function makes it possible to P(ause)/R(esume) or S(top) the continuer
%if npoint~=0 the Stop button is pressed en npoints becomes i

global sys driver_window calculation_progress;

npoints=0;
if sys.gui.pausespecial==1 %pause at special points
    if calculation_progress==2
%         set(status,'String','pause');
        waitfor(driver_window,'UserData',0);
        set(status,'String','computing');
        set(driver_window,'UserData',1);
    end
    if (index==(i-1)&&(i~=2))
%         set(status,'String','pause');
        figure(driver_window);
        set(driver_window,'UserData',1);
        waitfor(driver_window,'UserData',0);
        set(status,'String','computing');
    end
    if calculation_progress==0
        npoints = i;
    end      
elseif sys.gui.pauseeachpoint==1 %pause after each computed point
    if calculation_progress==0
        npoints = i;
    else
        set(status,'String','pause');
        set(driver_window,'UserData',1);
        waitfor(driver_window,'UserData',0);
        set(status,'String','computing');
    end
elseif sys.gui.pausenever==1 %pause only on command
    if calculation_progress==2
        set(status,'String','pause');
        set(driver_window,'UserData',1);
        waitfor(driver_window,'UserData',0);
        set(status,'String','computing');
    end
    if calculation_progress==0
        npoints = i;
    end      
end       
