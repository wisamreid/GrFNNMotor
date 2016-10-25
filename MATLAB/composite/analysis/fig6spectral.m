%% FIG6SPECTRAL
% This code plots Fig. 6 for the paper: spectra and 
% burst size.
% Version: 2013may23 by Cliff Kerr (cliffk@neurosim.downstate.edu)

%% Plotting
dosave=1;
figure('position',[100 100 400 600],'color','w')
plotparams={[95 77],[5 5],[15 -3],[1 1]};
kw={'wh','tc','bg','pa'};
fullkw={'WN','TC','BG','PD'};
linestyles={'--','--','-','-'};
popnames={'TC','IRE','E6','I6','IL6','E5B','E5R','I5','IL5','E4','I4','IL4','E2','I2','IL2'};
popinds=[   5,     6,  8,   9,  11,   12,   13,   14,  15,   16,  17,  18,   19,  21,  24];
for p=1:length(popnames), eval(sprintf('%s=%i;',popnames{p},popinds(p))); end

%% Plot labels
ax=axes('position',[0 0 1 1],'units','normalized'); hold on
text(0.03,0.96,'A','fontsize',16,'fontweight','bold')
text(0.03,0.48,'B','fontsize',16,'fontweight','bold')
xlim([0 1])
ylim([0 1])
set(ax,'visible','off')



%% Load data
nfiles=length(kw);
disp('Loading data...')
base='../output/output'; % Base file name -- needs to match runsim
locfiles=cell(nfiles,1);
lfpfiles=cell(nfiles,1);
spkfiles=cell(nfiles,1);
for f=1:nfiles
    fprintf('  File %i of %i...\n',f,nfiles)
    locfiles{f}=dlmread([base '-' kw{f} '-loc.txt']);
    lfpfiles{f}=dlmread([base '-' kw{f} '-lfp.txt']);
    spkfiles{f}=dlmread([base '-' kw{f} '-spk.txt']);
    lfpfiles{f}(:,1)=lfpfiles{f}(:,1)/1000; % Convert from ms to s
    spkfiles{f}(:,1)=spkfiles{f}(:,1)/1000; % Convert from ms to s
    spkfiles{f}(:,2)=spkfiles{f}(:,2)+1; % Change cell index to 1-based
    locfiles{f}(:,1)=locfiles{f}(:,1)+1; % Change cell index to 1-based
    ncells=length(locfiles{f}); % Number of cells in the model
end


%% Trim data
disp('Trimming data...')
trimends=2; % How many seconds of data to trim off the ends
for f=1:nfiles
    maxtime=lfpfiles{f}(end,1);
    start=trimends; 
    finish=maxtime-trimends;
    lfplimits=(lfpfiles{f}(:,1)>start & lfpfiles{f}(:,1)<finish);
    lfpfiles{f}=lfpfiles{f}(lfplimits,:);
    lfpfiles{f}(:,1)=lfpfiles{f}(:,1)-trimends;
    spklimits=(spkfiles{f}(:,1)>start & spkfiles{f}(:,1)<finish);
    spkfiles{f}=spkfiles{f}(spklimits,:);
    spkfiles{f}(:,1)=spkfiles{f}(:,1)-trimends;
end



%% Plot spectra
disp('Plotting spectra...')
cksubplot([1 2],[1 1],plotparams{:}); hold on
colors=vectocolor(1:nfiles+1)/1.3;
spectra=cell(nfiles,1);
layer=2;
rate=200;
smoothing=500;
maxfreq=100;
for f=1:nfiles
    thislfp=lfpfiles{f}(:,layer);
    thislfp=thislfp-mean(thislfp);
    thislfp=thislfp/mean(abs(thislfp));
    [spectra{f}, F]=spectrum('pl',thislfp,rate);
    spectra{f}(1)=spectra{f}(2); % To avoid problems with smoothing
    spectra{f}=cksmooth(spectra{f},smoothing);
end
for f=1:nfiles
    plot(F,spectra{f},'color',colors(f,:),'linewidth',2,'linestyle',linestyles{f})
end
set(gca,'yscale','log')
set(gca,'xscale','log')
axis tight
ylim([1e-3 1e2])
xlim([1 maxfreq])
legend(fullkw,'location','best')
xlabel('Frequency (Hz)')
ylabel('Power')



%% Plot burst frequency
disp('Plotting burst frequency...')

% Turn spike times into a raster
disp('  Binning spikes...')
maxtime=ceil(lfpfiles{1}(end,1));
binsize=0.01; % Bin size in s -- WARNING, results change based on the timestep!! :(
maxbins=ceil(maxtime/binsize)+1;
fullraster=zeros(nfiles,ncells,maxbins);
for f=1:nfiles
    thesespikes=round(spkfiles{f}(:,1)/binsize)+1; % Convert from time to bin index
    thesecells=spkfiles{f}(:,2);
    nspks=length(thesespikes);
    for s=1:nspks
        fullraster(f,thesecells(s),thesespikes(s))=fullraster(f,thesecells(s),thesespikes(s))+1;
    end
end
bursts=zeros(nfiles,maxbins);
bursthist=zeros(nfiles,ncells+1);
normalizedbursthist=zeros(nfiles,ncells+1);
burstsmooth=5;
count=0;
for f=1:nfiles
    fullraster(f,:,:)=fullraster(f,:,:)>1;
    bursts(f,:)=sum(fullraster(f,:,:),2);
    bursthist(f,:)=histc(bursts(f,:),0:ncells)/maxbins;
    bursthist(f,:)=cksmooth(bursthist(f,:),burstsmooth);
    bursthist(f,:)=bursthist(f,:)./sum(bursthist(f,:)); % This is to make sure that smoothing doesn't change the normalization
    spikeprob=sum(sum(fullraster(f,:,:))/numel(fullraster(f,:,:))); % The probability of a spike is just the number of spikes divided by the number of bins
    expectedhist=binopdf(0:ncells,ncells,spikeprob);
    normalizedbursthist(f,:)=bursthist(f,:)./expectedhist; % Express relative to shuffled data
end

% Plotting
cksubplot([1 2],[1 2],plotparams{:})
maxbursts=ncells;
for f=1:nfiles
    plot((0:maxbursts)/ncells*100,normalizedbursthist(f,:),'color',colors(f,:),'linewidth',2,'linestyle',linestyles{f}); hold on
end
set(gca,'yscale','log')
xlabel('Percent of cells firing simultaneously')
ylabel('Probability ratio')
xlims=[0 1.5];
xlim(xlims)
ylim([1e-2 1e75])
set(gca,'ytick',[1e0 1e25 1e50 1e75])
