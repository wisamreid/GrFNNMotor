function varargout = selectcycleConnecHom(varargin)

global path_sys gds HTHomds MC saddle_point driver_window dim_npos QU1; 
%dim_npos saddle_point QU1: from start_cont in gui_ConnecHom
prompt  = {'Enter the number of test intervals';'Enter the number of collocation points'};
title   = 'Choose ntst and ncol';
lines   = 1;
def     = {'40','4'};
answer  = inputdlg(prompt,title,lines,def);
if isempty(answer)||isempty(str2double(answer))
    ntst=40;
    ncol=4;
else
    ntst=str2double(answer{1});
    ncol=str2double(answer{2});    
end

filemat = strcat(gds.curve.new,'.mat');
filemat = fullfile(path_sys,gds.system,gds.diagram,filemat);
if exist(filemat,'file')
    load(filemat);
    if size(x,2)<4
        warndlg('No enough data available!');
        return;
    end
else
    warndlg('No enough data available!');
    return;
end    

ind = size(x,2); %number of points
while ((ind > 1) && isnan(x(end,ind))) || ((ind > 1) && (x(end,ind) - x(end,ind-1) >= 0))
    ind = ind-1;
end

%ind = 0: monotonous increasing
if ind <= 10
    %we take the whole orbit as starting orbit
    ind = size(x,2);
end
if ind == size(x,2)
    %monotonous decreasing
    verschil = x(end,2:end) - x(end,1:end-1);
    pos = find(verschil > 0);
    if size(pos,2) == 0
        warndlg('No enough data available!');
        return;
    end
end

x=x(1:gds.dim,1:ind);%omit eps1, just keep the coordinates
t=t(1:ind)';
tn=(t-t(1))/(t(end)-t(1));% tn(1) = 0 and tn(end)=1
str2='tempadwgyk_cycleconnechom.mat';
dir2 = fullfile(path_sys,gds.system,gds.diagram,str2);

[a,x,tn] = newmeshHT(x,tn,size(x,2)-1,1,ntst,ncol);

x=interp(tn,1,x,a,ncol);

HTHomds.msh = a;
s(1).data.timemesh = a;
s(1).data.ntst = ntst;
s(1).data.ncol = ncol;
s(1).data.parametervalues = param;
epsilon1 = norm(x(:,end)-saddle_point);

HTHomds.SParams = [];
for i = 1:dim_npos
    HTHomds.SParams{i,1} = strcat('SParam',num2str(i));    
    HTHomds.SParams{i,2} = 1/epsilon1*(x(:,end)-saddle_point)'*QU1(:,end-dim_npos+i);
end


HTHomds.P0 = param;
HTHomds.ntst = ntst;
HTHomds.ncol = ncol;
HTHomds.nphase = gds.dim;
x = reshape(x,size(x,2)*size(x,1),1);
HTHomds.T = (t(end)-t(1))/2;%half of the period
HTHomds.eps0 = gds.eps0;
HTHomds.eps1 = epsilon1;
s(1).data.T = HTHomds.T;
HTHomds.x0 = saddle_point;
HTHomds.extravec = [0 0 0];

HTHomds.UParams = [];    

if dim_npos >= 1
    HTHomds.UParams{1,1} = strcat('UParam',num2str(1));    
    HTHomds.UParams{1,2} = gds.c{1,2};        
end
if dim_npos >= 2
    HTHomds.UParams{2,1} = strcat('UParam',num2str(2));    
    HTHomds.UParams{2,2} = gds.c{2,2};    
end
if dim_npos >= 3
    for f = 3:dim_npos    
        HTHomds.UParams{f,1} = strcat('UParam',num2str(f));        
        HTHomds.UParams{f,2} = 0;        
    end
end

HTHomds.index = 0;
HTHomds.TestTolerance = gds.options.TestTolerance;
         
s(2)=[];
v=[];h=[];f=[];cds=[];ctype='ConnecHom ';point='ConnecHom ';num=1;
save(dir2,'x','v','s','h','f','cds','HTHomds','ctype','point','num')
    
    
index = 1;
gds.curve.old = 'tempadwgyk_cycleconnechom';feval(strcat('gui_',deblank(ctype)));
feval(gds.gui_load_point,index,x,'save',dir2,num);
s(num).label = ctype;
file = fullfile(path_sys,gds.system);
save(file,'gds');
list = get(MC.mainwindow.initial_point,'children');
tag = get(list,'Tag');
label = deblank(s(num).label);
for i = 1:length(list)
    [tag1,tag2] = strtok(tag{i},'_');
    if strcmp(label,tag1)
        type = strtok(tag2,'_');
        break
    else
        type = deblank(ctype);
    end
end
feval(strcat('gui_',type));
feval(gds.gui_point,tag{i});
load_matcont;
if size(MC.numeric_fig,1)==1; numeric;end