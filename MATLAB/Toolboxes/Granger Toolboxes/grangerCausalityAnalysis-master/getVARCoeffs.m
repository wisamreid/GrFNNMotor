function A = getVARCoeffs(bRaster, inds)


%% VAR
for i = 1:numel(inds)
    disp(strcat('Computing autoVAR for neuron # ',num2str(i), 'out of ',num2str(numel(inds))));
    results = vare(bRaster(inds(i),:),20);
    A{i} = results.nvar;
end