function Callback(varargin)
global gds path_sys;

h = gcbo;
tag = get(h,'Tag');
st  = str2double(strtok(tag,'edit'));
dim = gds.dim;
ndim = size(gds.parameters,1);
if (st==0)
    value = str2double(get(h,'String'));
    gds.time{1,2} = value;
    gds.integrator.tspan(2) = gds.integrator.tspan(2)-gds.integrator.tspan(1)+gds.time{1,2};
    gds.integrator.tspan(1) = gds.time{1,2};
    file = fullfile(path_sys,gds.system);save(file,'gds');    
    return
end
if (st <= dim) % case coordinates
     value = str2double(get(h,'String'));
     gds.coordinates{st,2} = value;
elseif ((st <= dim+ndim)&&(st > dim)) %case parameters
    i = st-gds.dim;
    value = str2double(get(h,'String'));
    gds.parameters{i,2}=value;
elseif ((st <= 2*dim+ndim)&&(st>dim+ndim))
     i = st-(dim+ndim);
     value = str2double(get(h,'String'));
     gds.x0{i,2} = value;     
elseif ((st <= 3*dim+ndim)&&(st>2*dim+ndim)) 
     i = st-(2*dim+ndim);
     value = str2double(get(h,'String'));
     gds.x1{i,2} = value;
elseif strcmp(tag,'cc1')
    value = str2double(get(h,'String'));
    gds.c{1,2}=value;
elseif strcmp(tag,'cc2')
    value = str2double(get(h,'String'));    
    gds.c{2,2}=value;            
elseif strcmp(tag,'eps0val')
    value = str2double(get(h,'String'));    
    gds.eps0=value;   
else
    gds.period = str2double(get(h,'String'));
end
file = fullfile(path_sys,gds.system);
save(file,'gds');  