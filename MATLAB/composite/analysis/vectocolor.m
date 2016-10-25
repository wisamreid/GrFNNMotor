%% VECTOCOLOR
% This function converts a vector of N values into an Nx3 matrix of color
% values according to the current colormap. It automatically scales the 
% vector to provide maximum dynamic range for the color map.
%
% Usage:
%   colors=vectocolor(vector)
%
% where:
%   colors is an Nx3 vector of RGB color values
%   vector is the input vector
%
% Example:
%	n=1000;
%	x=randn(n,1);
%	y=randn(n,1);
%	z=log(x.^2+y.^2);
%	c=vectocolor(z);
%	scatter3(x,y,z,20,c,'filled')
%
% Version: 2012jul03 by Cliff Kerr (cliffk@neurosim.downstate.edu)

function colors=vectocolor(vector)

if isempty(get(0,'CurrentFigure')), tmpfig=figure; end % If no figure exists, create one
map=colormap; % Get the current color map
vector=round(normalize(vector,[1 length(map)])); % Convert to integers for the colormap
colors=map(vector,:); % Pick colors for each element in the vector
if exist('tmpfig','var'); close(tmpfig); end; % If created figure, close it

end





%% NORMALIZE
% This function takes a matrix of any size and changes its range to the
% specified value.
% 
% Usage:
%     newmatrix=normalize(oldmatrix,limits)
% 
% where:
%     newmatrix is the normalized matrix
%     oldmatrix is the original matrix
%     limits is a 2-element vector of upper and lower limits (optional; default [0,1])
% 
% Example:
%     oldvector=randn(1000,1);
%     limits=[-1 1];
%     newvector=normalize(oldvector,limits);
%     hist(newvector,100)
%     
% Version: 2012oct03 by cliffk

function newmatrix=normalize(varargin)

% Handle input arguments
oldmatrix=varargin{1}; % Matrix to be normalized
if nargin>1, limits=varargin{2}; else limits=[0,1]; end % Range to normalize matrix to
if numel(oldmatrix)==0 % Empty array, so return empty
    newmatrix=[];
elseif numel(oldmatrix)==1 % Single-element array, so return midpoint of limits
    newmatrix=sum(limits)/2;
else % All other cases, so do proper computation
    newmatrix=oldmatrix; % Copy matrix
    newmatrix=newmatrix-min(newmatrix(:)); % Set minimum to 0
    newmatrix=newmatrix/max(newmatrix(:)); % Set maximum to 1
    newmatrix=newmatrix*(range(limits)); % Get range right
    newmatrix=newmatrix+limits(1); % Get minimum right
end

end
