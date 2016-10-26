function status = integ_prs(t,y,flag,varargin)
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

  
global MC driver_window calculation_progress;
status = 0;
if isempty(flag) % odephas2(t,y) [v5 syntax] or odephas2(t,y,'')
    if calculation_progress==0
       status = 1;
    else 
       status = 0;
    end;
    if calculation_progress==2
        set(MC.mainwindow.mstatus,'String','pause');
        waitfor(driver_window,'UserData',[0]);
        set(MC.mainwindow.mstatus,'String','computing');
    end;       
else
    switch(flag)
    case 'init'                           % odephas2(tspan,y0,'init')
        set(MC.mainwindow.mstatus,'String','computing');
        stop;
    case 'done'         % integ_plot([],[],'done')   
        set(MC.mainwindow.mstatus,'String','ready');
        if ishandle(driver_window),delete(driver_window);end
    end
end
drawnow;