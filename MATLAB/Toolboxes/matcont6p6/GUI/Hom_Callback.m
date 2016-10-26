function Hom_Callback(varargin)
global gds path_sys;
h = gcbo;
tag = get(h,'Tag');
st  = str2double(strtok(tag,'edit'));
value = str2double(get(h,'String'));
if st == 1 %case parameters
    gds.T=value;
elseif st == 2
    gds.eps0 = value;
else
    gds.eps1 = value;
end
file = fullfile(path_sys,gds.system);
save(file,'gds');  

