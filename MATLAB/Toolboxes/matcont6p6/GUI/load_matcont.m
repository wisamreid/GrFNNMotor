function load_matcont(varargin)
%function sets the point type etc. in the main(matcont)-window
global gds MC calculation_progresss
for i=1:5
    for j=1:4
        if (gds.der(j,i)==1)
            switch j
            case 1
                str(i)='N';
            case 3
                str(i)='R';
            case 4
                if (exist('sym')==2)
                    str(i)='S';
                else
                    str(i)='F';
                end
            end
        end
    end
end
set(MC.mainwindow.msystem,'String',gds.system)
set(MC.mainwindow.mcurve,'String',gds.curve.new)
set(MC.mainwindow.mpointtype,'String',gds.point)
set(MC.mainwindow.mcurvetype,'String',gds.type)
set(MC.mainwindow.mderivatives,'String',str);
set(MC.mainwindow.duration,'string','');
set(MC.mainwindow.mstatus,'string','');

