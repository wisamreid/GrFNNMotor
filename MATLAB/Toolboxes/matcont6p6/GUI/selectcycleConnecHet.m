function varargout = selectcycleConnecHet(varargin)

global path_sys gds HTHetds MC driver_window dim_npos dim_nneg QS; 

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

ind = size(x,2); %aantal punten
while (ind > 1) && (x(end,ind) - x(end,ind-1) >= 0)
    ind = ind-1;
end

%ind = 0: monotoon stijgend
if ind <= 10
    ind = size(x,2);
end
if ind == size(x,2)
   %monotoon dalend
    verschil = x(end,2:end) - x(end,1:end-1);
    pos = find(verschil > 0);
    if size(pos,2) == 0
        warndlg('No enough data available!');
        return;
    end
end

x=x(1:gds.dim,1:ind);%hier laat je de eps1 weg, hou je enkel de coordinaten over
t=t(1:ind)';
tn=(t-t(1))/(t(end)-t(1));%zo is tn(1) = 0 en tn(end)=1
str2='tempadwgyk_cycleconnechet.mat';
dir2 = fullfile(path_sys,gds.system,gds.diagram,str2);
[a,x,tn] = newmeshHT(x,tn,size(x,2)-1,1,ntst,ncol);
x = interp(tn,1,x,a,ncol);

s(1).data.timemesh = a;
s(1).data.ntst = ntst;
s(1).data.ncol = ncol;
s(1).data.parametervalues = param;
epsilon1 = norm(x(:,end)-vertcat(gds.x1{:,2}));

HTHetds.SParams = [];
HTHetds.UParams = [];
for i = 1:(gds.dim-dim_nneg)
    HTHetds.SParams{i,1} = strcat('SParam',num2str(i));    
    HTHetds.SParams{i,2} = 1/epsilon1*(x(:,end)-vertcat(gds.x1{:,2}))'*QS(:,dim_nneg+i);
end

HTHetds.P0 = param;
HTHetds.ntst = ntst;
HTHetds.ncol = ncol;
HTHetds.nphase = gds.dim;
x = reshape(x,size(x,2)*size(x,1),1);%je maakt er hier 1 lange vector van
HTHetds.T = (t(end)-t(1))/2;%de helft van de periode
HTHetds.eps0 = gds.eps0; %hier voeg je epsilon0 toe
HTHetds.eps1 = epsilon1;
s(1).data.T = HTHetds.T;
HTHetds.x0 = vertcat(gds.x0{:,2});
HTHetds.x1 = vertcat(gds.x1{:,2});
HTHetds.extravec = [0 0 0];
    
if dim_npos >= 1
    HTHetds.UParams{1,1} = strcat('UParam',num2str(1));    
    HTHetds.UParams{1,2} = gds.c{1,2};        
end
if dim_npos >= 2
    HTHetds.UParams{2,1} = strcat('UParam',num2str(2));    
    HTHetds.UParams{2,2} = gds.c{2,2};    
end
if dim_npos >= 3
    for f = 3:dim_npos    
        HTHetds.UParams{f,1} = strcat('UParam',num2str(f));        
        HTHetds.UParams{f,2} = 0;        
    end
end

HTHetds.index = 0;
HTHetds.TestTolerance = gds.options.TestTolerance;
         
s(2)=[];
v=[];h=[];f=[];cds=[];ctype='ConnecHet ';point='ConnecHet ';num=1;
save(dir2,'x','v','s','h','f','cds','HTHetds','ctype','point','num')
    
    
index = 1;
gds.curve.old = 'tempadwgyk_cycleconnechet';feval(strcat('gui_',deblank(ctype)));
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