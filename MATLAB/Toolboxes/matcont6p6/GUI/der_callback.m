function der_callback
global gds;
h = gcbo;
num = get(h,'UserData');
str = get(h,'Tag');
if (num==1)
    if strcmp(str,'editp')
     gds.jacp = get(h,'String');
    else 
     gds.jac = get(h,'String');
    end
end
if (num==2)
    if strcmp(str,'editp')
     gds.hessp = get(h,'String');
    else 
     gds.hess = get(h,'String');
    end
end

