% [bRaster,mask] = makeBurstRaster(asdf)
%
% Author : Madhavun Candadai Vasu
% Date   : Mar 10, 2016
%
% This function does the following:
%   - intersect mask with asdf_raw to create burst asdf
%   - convolve raster to convert binary raster to continuous signal


function bRaster = makeBurstRaster(raster,mask)

% find start and end of each burst from mask
fst = []; 
for i=1:numel(mask)-1
    fst(end+1) = mask(i+1)-mask(i);
end
disp('Generated masks');
maskStarts = find(fst == 1);
maskEnds = find(fst == -1);
disp(strcat('Number of masks - ',num2str(numel(maskStarts))));

% check for size match
if numel(maskStarts) ~= numel(maskEnds)
    error('Error message: number of mask startings points is not equal to number of mask ending points');
    return
end

% create empty bAsdf
for i=1:N
    bAsdf{i} = [];
end
bRaster = [];
mt = 0; %max Time unit
% convert each mask area to raster and generate bAsdf
for i=1:numel(maskStarts)
    clc;
    disp('Generated masks');
    disp(strcat('Number of masks - ',num2str(numel(maskStarts))));
    disp(strcat('Progress - ',num2str(i/numel(maskStarts)*100),'%'));
    cutasdf = ASDFChooseTime(asdf,maskStarts(i),maskEnds(i));
    [raster, ~] = ASDFToSparse(cutasdf);
    raster = full(raster);
%     % generate bAsdf
%     if mt ~= 0
%         mt = findMaxTimeUnit(bAsdf);
%     end
%     for n=1:N
%         bAsdf{i} = [bAsdf{i},mt+find(raster(n,:) == 1)]
%     end
    bRaster = [bRaster raster];
end
% bAsdf{N+1} = 1;
% bAsdf{N+2} = findMaxTimeUnit(bAsdf);

disp('Done!');