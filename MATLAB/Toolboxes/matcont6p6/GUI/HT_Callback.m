function HT_Callback(varargin)
global gds path_sys

h = gcbo;
tag = get(h,'Tag');
st  = str2double(strtok(tag,'edit'));

value = str2double(get(h,'String'));
if ((st <= gds.dim+size(gds.parameters,1)+size(gds.UParams,1)) && (st>gds.dim+size(gds.parameters,1)))%case UParams (c's) bij connection
    i = st-gds.dim-size(gds.parameters,1);
    value = str2double(get(h,'String'));
    gds.UParams{i,2}=value;
elseif ((st <= gds.dim+size(gds.parameters,1)+size(gds.UParams,1)+size(gds.SParams,1)) && (st>gds.dim+size(gds.parameters,1)+size(gds.UParams,1)))%case SParams (c's) bij connection
    i = st-gds.dim-size(gds.parameters,1)-size(gds.UParams,1);
    value = str2double(get(h,'String'));
    gds.SParams{i,2}=value;
elseif st == 999
    gds.T = value;
elseif st == 10000
    gds.eps1 = value;
elseif st == 1001
    gds.eps1tol = value;
end

file = fullfile(path_sys,gds.system);
save(file,'gds');  