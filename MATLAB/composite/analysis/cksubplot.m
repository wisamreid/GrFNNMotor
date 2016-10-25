%% CKSUBPLOT
% This function is a generalization of subplot that makes it easier to
% change the space between panels.
%
% Usage:
%   axishandle=cksubplot([cols,rows],[thiscol,thisrow],[horizsize%,vertsize%],[xmargin%,ymargin%],[xoffset%,yoffset%],[xlabels,ylabels])
% 
% Example:
%  for i=1:5, for j=1:2, h=cksubplot([5,2], [i,j], [90,90], [5,5], [10,10], [0,0]), end, end
%
% Version: 2012aug16 by Cliff Kerr (cliffk@neurosim.downstate.edu)

function varargout=cksubplot(panellayout,thispanel,varargin)

% Set defaults: completely fill the figure
if nargin>=3, panelsize=varargin{1}; else panelsize=[100,100]; end
if nargin>=4, figmargin=varargin{2}; else figmargin=[0,0]; end
if nargin>=5, figoffset=varargin{3}; else figoffset=[0,0]; end
if nargin>=6, xylabels=varargin{4}; else xylabels=[0,0]; end
if nargin>=7, otherargs=varargin{5:end}; else otherargs={}; end

% Handle scalar thispanel input
if numel(thispanel)==1
    col=mod(thispanel-1,panellayout(1))+1;
    row=floor((thispanel-1)/panellayout(1))+1;
    thispanel=[col row];
end

% Convert from percentages to decimals
panelsize=panelsize/100;
figmargin=figmargin/100;
figoffset=figoffset/100;

% Convert user entries into [x,y] coordinates
width  =panelsize(1)*(1-figmargin(1))*(1-figoffset(1))/panellayout(1);
height =panelsize(2)*(1-figmargin(2))*(1-figoffset(2))/panellayout(2);
xoffset=figoffset(1)+(thispanel(1)-1)*(1-figmargin(1))/panellayout(1);
yoffset=1-(figoffset(2)+(thispanel(2))*(1-figmargin(2))/panellayout(2));

% Create panel
axishandle=axes('units','normalized','position',[xoffset yoffset width height],otherargs{:}); % Add any extra inputs here
if xylabels(1)==0 && thispanel(2)~=panellayout(2), set(axishandle,'xticklabel',[]), end % Turn off interior x-axis labels
if xylabels(2)==0 && thispanel(1)~=1, set(axishandle,'yticklabel',[]), end % Turn off interior y-axis labels
box on % I almost always want this, so make it the default


if nargout>0, varargout{1}=axishandle; end


end
