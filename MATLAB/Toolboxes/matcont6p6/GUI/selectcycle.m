function varargout = selectcycle(varargin)

global path_sys gds MC driver_window;
prompt  = {'Enter the tolerance to select a cycle';'Enter the number of test intervals'};
title   = 'Choose tolerance and ntst';
lines   = 1;
def     = {'1e-2','20'};
answer  = inputdlg(prompt,title,lines,def);
if isempty(answer)||isempty(str2double(answer{1}))
    tolerance = 1e-2;
    ntst=20;
else
    tolerance = str2double(answer{1});
    ntst=str2double(answer{2});
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
xstart=round(size(x,2)/3);
amin=sum(abs(x(:,xstart:end)-x(:,1)*ones(1,size(x(:,xstart:end),2))));
ep=tolerance;
amin(find(amin(:)<ep))=inf;
[pp,qq]=min(amin-(1e-4));
if (pp==Inf)|pp<0
   ind=[]; 
else
    ind = qq(1)+xstart-1;
end
if isempty(ind)
    warndlg('No cycle can be found!')
    return; 
end
x=x(:,1:ind);
t=t(1:ind)';
tn=(t-t(1))/(t(end)-t(1));
str2='tempadwgyk_cycle.mat';
dir2 = fullfile(path_sys,gds.system,gds.diagram,str2);
[a,x,tn] = newmeshcycle(x,tn,size(x,2)-1,1,ntst,4);
x = interp(tn,1,x,a,4);
s(1).data.timemesh = a;
s(1).data.ntst = ntst;
s(1).data.ncol = 4;
s(1).data.parametervalues = param;
lds.P0 = param;
lds.ActiveParams = 1;
lds.ntst = ntst;
lds.ncol = 4;
x = reshape(x,size(x,2)*size(x,1),1);
x(end+1) = t(end)-t(1);
s(1).data.T = x(end);
lds.T = x(end);
x(end+1) = param(1);
s(2)=[];
v=[];h=[];f=[];cds=[];ctype='LC ';point='LC ';num=1;
save(dir2,'x','v','s','h','f','cds','lds','ctype','point','num')
index = 1;
gds.curve.old = 'tempadwgyk_cycle';feval(strcat('gui_',deblank(ctype)));
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