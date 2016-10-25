%% CKSMOOTH
% This function applies a smoothing kernel to a 1- or 2-dimensional array.
% Be warned that it is an incredibly inefficient implementation!
%
% Usage:
%   B=cksmooth(A,[n],[dim])
%
% where
%   B is a smoothed vector/array with the same dimensions as the original
%   A is the original vector or array
%   n is the number of adjacent elements to smooth over; default 5
%   dim is the dimension to sum over (for 2D arrays)
%
% Version: 2012feb18 by Cliff Kerr (cliffk@neurosim.downstate.edu)

function output=cksmooth(varargin)

% Handle input
repeats=5; % How many repeats by default?
dim=0; % Which dimensions to smooth over?
if nargin>=1, input=varargin{1}; end % First argument is input array
if nargin>=2 && ~isempty(varargin{2}), repeats=varargin{2}; end
if nargin>=3, dim=varargin{3}; end
if nargin<1 || nargin>=4, error('Function takes 1-3 input arguments.'); end
origsize=size(input);
if ndims(input)>2, error('You want to smooth a 3D array!? I don''t think so, buddy!'), end

% Define kernels
kernel1=[1 2 1]; % Define a 1D smoothing kernel
kernel2=[1 2 1; 2 4 2; 1 2 1]; % Column/row 2D smoothing kernel
kernel1=kernel1/sum(kernel1(:)); % Normalize
kernel2=kernel2/sum(kernel2(:)); % Normalize

% 1D smoothing
if numel(input)==length(input) % It's a vector
    output=reshape(input,length(input),1); % Turn into a column vector
    for i=1:repeats
        padded=[output(1+1);output;output(end-1)]; % Reflect the edges
        output=conv(padded,kernel1,'valid'); % Perform 1D convolution
    end

% 2D smoothing
else % It's a matrix
    output=input;
    if dim==0 % Smooth over both dimensions
        padded=zeros(size(output,1)+2,size(output,2)+2); % Make a padded array
        for i=1:repeats
            padded(2:end-1,2:end-1)=output; % Put array in the middle
            padded([1 end],2:end-1)=padded([2 end-1],2:end-1); % Add the top and bottom
            padded(:,[1 end])=padded(:,[2 end-1]); % Add the sides
            output=conv2(padded,kernel2,'valid'); % Perform 2D convolution
        end
    elseif dim==1 % Smooth over columns
        padded=zeros(size(output,1)+2,size(output,2)); % Make a padded array
        for i=1:repeats
            padded(2:end-1,:)=output; % Put array in the middle
            padded([1 end],:)=padded([2 end-1],:); % Add the top and bottom
            output=conv2(padded,kernel1','valid'); % Perform 2D convolution using the transposed 1D kernel
        end
    elseif dim==2 % Smooth over rows
        padded=zeros(size(output,1),size(output,2)+2); % Make a padded array
        for i=1:repeats
            padded(:,2:end-1)=output; % Put array in the middle
            padded(:,[1 end])=padded(:,[2 end-1]); % Add the sides
            output=conv2(padded,kernel1,'valid'); % Perform 2D convolution using the 1D kernel
        end
    else
        error('Invalid dimension to smooth over: must be blank (0), 1, or 2')
    end
end

output=reshape(output,origsize); % Make sure it has the original dimensions

end
