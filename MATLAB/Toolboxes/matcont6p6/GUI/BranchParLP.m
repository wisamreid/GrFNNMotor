function BranchParLP(varargin)
% Branch_Par Application M-file 
% last change on 20/08/03 15:30  
global gds path_sys MC
h = gcbo;
tag  = get(h,'Tag');
file = fullfile(path_sys,gds.system);
m = get(h,'Value');
user = get(h,'UserData');num=user.num;
if ~isfield(gds,'BranchParams'),gds.BranchParams=[];end
if isempty(gds.BranchParams)
    st  = [];
else
    st  = gds.BranchParams;
end
switch m
    case 0 % yes
        x = find(st==str2num(num));
        if isempty(x)
            return
        else
            st(x) = '';str = mat2str(st);
            val = sort(str2num(str2mat(str)));
            gds.BranchParams = val;
            save(file,'gds');
        end            
    case 1 % no        
        if isempty(gds.BranchParams)
            str = '[]';
        else
            str = mat2str(gds.BranchParams);
        end
        if (size(gds.BranchParams,2)==1)
            val = sort(str2num(str2mat(sprintf('[%s%c%s]',str,char(32),num))));
            gds.BranchParams = val;
            save(file,'gds');
        else
            str(end) = '';
            val = sort(str2num(str2mat(sprintf('%s%c%s]',str,char(32),num))));
            gds.BranchParams = val;
            save(file,'gds');
        end
end
feval(gds.gui_singularities);
if ~isempty(MC.numeric_fig),numeric;end
