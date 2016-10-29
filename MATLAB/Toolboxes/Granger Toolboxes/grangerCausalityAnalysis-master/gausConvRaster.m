% z = gausConvRaster(asdf,startT,endT)
%
% Author : Madhavun Candadai Vasu
% Date   : Feb 25, 2016
%
% This function does the following:
%   - takes in raw asdf of entire recording
%   - cuts it according to specified start and end times
%   - converts to raster format 
%   - convolves every neuron's raster with a gaussian
%   - returns convolved raster, z

function z = gausConvRaster(fraster,startT,endT)

nNeurons = size(fraster,1);

%% prep data
% cutasdf = ASDFChooseTime(asdf,startT,endT);
% [raster,binunit] = ASDFToSparse(cutasdf);
% fraster = full(raster);
%image(fraster*100);

%% create Gaussian kernel
sigma=1;
x=[-sigma*10:sigma/50:sigma*10];
gk = exp(-(x/sigma).^2/2)/(sigma*sqrt(2*pi));
gk = gk.*100; % gain

%% convolve
for neuronId = 1:nNeurons
    %choose data
    d = fraster(neuronId,:);
    
    if sum(d) > 0
        % convolve and %plot
        dgk = conv(d,gk);
        z(neuronId,:) = dgk(ceil(length(gk)/2):end-floor(length(gk)/2)); % Alignment
    else
        z(neuronId,:) = zeros(1,size(fraster,2));
    end
end

%% normalize
for i=1:size(z,1)
    if sum(z(i,:)) > 0
        z(i,:) = z(i,:)./max(z(i,:));
    end
end

%% plot
%figure; plot(gk);
%
% n = 2;
% figure; hold on;
% stem(fraster(n,:))
% plot(z(n,:),'k')
% title(strcat('Neuron index - ',num2str(n)))
% 
% n = 25;
% figure; hold on;
% stem(fraster(n,:))
% plot(z(n,:),'k')
% title(strcat('Neuron index - ',num2str(n)))
% 
% n = 14;
% figure; hold on;
% stem(fraster(n,:))
% plot(z(n,:),'k')
%title(strcat('Neuron index - ',num2str(n)))
