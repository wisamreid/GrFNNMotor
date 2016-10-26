function varargout = SBplot(varargin)
% SBplot - plots given data.
%
% USAGE:
% ======
% [] = SBplot(t,values)
% [] = SBplot(t,values,names)
% [] = SBplot(t,values,names,legendtext)
% [] = SBplot(t,values,names,legendtext,marker)
%
% t: column vector with time information
% values: matrix with data where each row corresponds to one time point and
%    each column to a different variable
% names: cell-array with the names of the data variables
% legendtext: cell-array of same length as names with text to be used for
%   the legend.
% marker: marker and line style for plot
%
% DEFAULT VALUES:
% ===============
% names: the plotted variables obtain the name 'x1', 'x2', ...
% legendtext: same as names
% marker: '-'

% Information:
% ============
% Systems Biology Toolbox for MATLAB
% Copyright (C) 2005-2007 Fraunhofer-Chalmers Research Centre
% Author of the toolbox: Henning Schmidt, henning@hschmidt.de
% 
% This program is free software; you can redistribute it and/or
% modify it under the terms of the GNU General Public License
% as published by the Free Software Foundation; either version 2
% of the License, or (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details. 
% 
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SBplot_OpeningFcn, ...
                   'gui_OutputFcn',  @SBplot_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end
if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT

% --- Executes just before SBplot is made visible.
function SBplot_OpeningFcn(hObject, eventdata, handles, varargin)
% determine number of arguments passed to SBplot by the user
nvarargin = nargin-3;
% process variable arguments
if nvarargin == 2,
    time = varargin{1};
    data = varargin{2};
    dataNames = {};
    for k = 1:size(data,2);
        dataNames{k} = sprintf('x%d',k);
    end
    legendtext =  dataNames;
    marker = '-';
elseif nvarargin == 3,
    time = varargin{1};
    data = varargin{2};
    dataNames = varargin{3};
    legendtext =  dataNames;
    marker = '-';
elseif nvarargin == 4,
    time = varargin{1};
    data = varargin{2};
    dataNames = varargin{3};
    legendtext =  varargin{4};
    marker = '-';
elseif nvarargin == 5,
    time = varargin{1};
    data = varargin{2};
    dataNames = varargin{3};
    legendtext =  varargin{4};
    marker = varargin{5};
else 
    error('Wrong number of input arguments');
end
% Check data consistency
if size(time,1) ~= size(data,1),
    error('Different number of time points and time points in data.');
end
if length(dataNames) ~= size(data,2),
    error('Different number of variable data and variable names.');
end
set(handles.xaxisselection,'String',{'TIME',dataNames{:}});
set(handles.yaxisselection,'String',dataNames);
% Set all the plot data also in the handles structure to be accessed by 
% all callback functions
handles.time = time;
handles.data = data;
handles.dataNames = dataNames;
handles.legentext = legendtext;
handles.marker = marker;
handles.dataPlotType = 'plot';
% Initialize export figure handle
handles.exportFigureHandle = [];
handles.grid = 0;
% Doing a first plot
doPlot(handles);
% Choose default command line output for SBplot
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
return

% --- Outputs from this function are returned to the command line.
function varargout = SBplot_OutputFcn(hObject, eventdata, handles) 
varargout{1} = handles.output;
return

% --- Executes on button press in zoombutton.
function zoombutton_Callback(hObject, eventdata, handles)
% toogle the zoom in the figure
zoom
return

% --- Executes on selection change in xaxisselection.
function xaxisselection_Callback(hObject, eventdata, handles)
try
    doPlot(handles);
catch
    errordlg('This selection is not possible.','Error','on');               
end
return

% --- Executes on selection change in yaxisselection.
function yaxisselection_Callback(hObject, eventdata, handles)
try
    doPlot(handles);
catch
    errordlg('This selection is not possible.','Error','on');               
end
return

% --- Executes on button press in gridbutton.
function gridbutton_Callback(hObject, eventdata, handles)
% toogle the grid in the figure
grid
if handles.grid == 1,
    handles.grid = 0;
else
    handles.grid = 1;
end
% Update handles structure
guidata(hObject, handles);
return

% --- Executes on button press in plot.
function plot_Callback(hObject, eventdata, handles)
handles.dataPlotType = 'plot';
% Update handles structure
guidata(hObject, handles);
doPlot(handles);
return

% --- Executes on button press in loglog.
function semilogx_Callback(hObject, eventdata, handles)
handles.dataPlotType = 'semilogx';
% Update handles structure
guidata(hObject, handles);
doPlot(handles);
return

% --- Executes on button press in semilogx.
function semilogy_Callback(hObject, eventdata, handles)
handles.dataPlotType = 'semilogy';
% Update handles structure
guidata(hObject, handles);
doPlot(handles);
return

% --- Executes on button press in loglog.
function loglog_Callback(hObject, eventdata, handles)
handles.dataPlotType = 'loglog';
% Update handles structure
guidata(hObject, handles);
doPlot(handles);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXPORT FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function export_Callback(hObject, eventdata, handles)
if isempty(handles.exportFigureHandle),
    figH = figure;
    handles.exportFigureHandle = figH;
    % Update handles structure
    guidata(hObject, handles);
else
    figH = handles.exportFigureHandle;
    figure(figH);
end
nrow = str2num(get(handles.nrow,'String'));
ncol = str2num(get(handles.ncol,'String'));
nnumber = str2num(get(handles.nnumber,'String'));
subplot(nrow,ncol,nnumber);
doPlot(handles);
if handles.grid == 1,
    grid;
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REQUEST NEW EXPORT FIGURE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function newexportfigure_Callback(hObject, eventdata, handles)
handles.exportFigureHandle = [];
% Update handles structure
guidata(hObject, handles);
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PLOT FUNCTION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function doPlot(handles)
time = handles.time;
data = handles.data;
dataNames = handles.dataNames;
xaxis = handles.xaxisselection;
yaxis = handles.yaxisselection;
% get variable that is chosen for the x-axis
indexX = get(xaxis,'Value');
% get variables that are chosen for the y-axis
indexY = get(yaxis,'Value');
if indexX == 1,
    if size(time,2) == 1,
        xvariable = time;
    else
        xvariable = time(:,indexY);
    end
else
    xvariable = data(:,indexX-1);
end
yvariables = data(:,indexY);
yvariablesNames = dataNames(indexY);
feval(handles.dataPlotType,xvariable,yvariables,handles.marker)
legend(handles.legentext(indexY)); 
% write axis labels
if indexX == 1,
    xlabel('Time');
else
    xlabel(dataNames(indexX-1));
end
return



