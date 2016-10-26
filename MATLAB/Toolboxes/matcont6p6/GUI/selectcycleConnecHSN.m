 function varargout = selectcycleConnecHSN(varargin)

global path_sys gds HTHSNds MC driver_window dim_npos QS1 eig0; 
%dim_npos saddle_point QU1: komt van start_cont in gui_Con
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

%ind = size(x,2);

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
str2='tempadwgyk_cycleconnechsn.mat';
dir2 = fullfile(path_sys,gds.system,gds.diagram,str2);
[a,x,tn] = newmeshHT(x,tn,size(x,2)-1,1,ntst,ncol);
x = interp(tn,1,x,a,ncol);

s(1).data.timemesh = a;
s(1).data.ntst = ntst;
s(1).data.ncol = ncol;
s(1).data.parametervalues = param;
epsilon1 = norm(x(:,end)-vertcat(gds.x0{:,2}));
%vanaf hier moet je een onderscheid maken

HTHSNds.SParams = [];
HTHSNds.UParams = [];
dim_nneg = gds.dim-dim_npos-1;
Qtot = [QS1(:,1:dim_nneg) eig0];
[Qn,Rn] = qr(Qtot);
for i = 1:dim_npos
    HTHSNds.SParams{i,1} = strcat('SParam',num2str(i));    
    HTHSNds.SParams{i,2} = 1/epsilon1*(x(:,end)-vertcat(gds.x0{:,2}))'*Qn(:,dim_nneg+1+i);
end


HTHSNds.P0 = param;
HTHSNds.ntst = ntst;
HTHSNds.ncol = ncol;
HTHSNds.nphase = gds.dim;
x = reshape(x,size(x,2)*size(x,1),1);%je maakt er hier 1 lange vector van
HTHSNds.T = (t(end)-t(1))/2;%de helft van de periode
HTHSNds.eps0 = gds.eps0; %hier voeg je epsilon0 toe
HTHSNds.eps1 = epsilon1;
s(1).data.T = HTHSNds.T;
HTHSNds.x0 = vertcat(gds.x0{:,2});
HTHSNds.extravec = [0 0 0];

HTHSNds.UParams{1,1} = strcat('UParam',num2str(1));    
HTHSNds.UParams{1,2} = gds.c{1,2};        

if dim_npos >= 1
    for f = 2:dim_npos+1    
        HTHSNds.UParams{f,1} = strcat('UParam',num2str(f));            
        HTHSNds.UParams{f,2} = 0;            
    end
end

HTHSNds.index = 0;
HTHSNds.TestTolerance = gds.options.TestTolerance;
         
s(2)=[];
v=[];h=[];f=[];cds=[];ctype='ConnecHSN ';point='ConnecHSN ';num=1;
save(dir2,'x','v','s','h','f','cds','HTHSNds','ctype','point','num')
       
index = 1;
gds.curve.old = 'tempadwgyk_cycleconnechsn';feval(strcat('gui_',deblank(ctype)));

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