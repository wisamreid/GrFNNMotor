%% FIG5FIRING
% This code plots Fig. 5 for the paper: a raster, time series, average
% firing rates, and Fano factors.
% Version: 2013may23 by Cliff Kerr (cliffk@neurosim.downstate.edu)

%% Define things
rng(1)
figure('position',[100 100 800 800],'color','w')
popnames={'TCR','TRN','E6','I6','IL6','E5B','E5R','I5','IL5','E4','I4','IL4','E2','I2','IL2'};
popinds=[   5,     6,  8,   9,  11,   12,   13,   14,  15,   16,  17,  18,   19,  21,  24];
for p=1:length(popnames), eval(sprintf('%s=%i;',popnames{p},popinds(p))); end
popinds2=zeros(max(popinds),1);
for p=1:length(popinds)
    popinds2(popinds(p))=p;
end
namevec=cell(max(popinds),1);
for p=1:length(popinds), namevec{popinds(p)}=popnames{p}; end
popcolors=flipud([[0.42,0.67,0.84];[0.50,0.80,1.00];[0.90,0.32,0.00];[0.34,0.67,0.67];[0.42,0.82,0.83];[0.90,0.59,0.00];[0.33,0.67,0.47];[0.42,0.83,0.59];[0.90,0.76,0.00];[1.00,0.85,0.00];[0.71,0.82,0.41];[0.57,0.67,0.33];[1.00,0.38,0.60];[0.5,0.2,0.0];[0.0,0.2,0.5]]); % Colors for each cell population


%% Plot labels
ax=axes('position',[0 0 1 1],'units','normalized'); hold on
text(0.03,0.97,'A','fontsize',16,'fontweight','bold')
text(0.56,0.97,'B','fontsize',16,'fontweight','bold')
text(0.03,0.28,'C','fontsize',16,'fontweight','bold')
text(0.55,0.28,'D','fontsize',16,'fontweight','bold')
xlim([0 1])
ylim([0 1])
set(ax,'visible','off')


%% Load data
kw={'wh','tc','bg','pa'};
fullkw={'WN','TC','BG','PD'};

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
trimends=3; % How many seconds of data to trim off the ends
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


%% Plot raster
cksubplot([2 2], [1 1], [80,130], [0,0], [11,12], [0,0])
thismodel=3; % Which model to use: 3 = healthy basal ganglia
xmax=5;
fractoplot=0.2; % What fraction of spikes toplot
labeloffset=10;
TCRoffset=15;
cellids=locfiles{thismodel}(:,1);
cellpops=locfiles{thismodel}(:,2);
plotoffset=0;
for j=length(popinds):-1:1
    alltime=spkfiles{thismodel}(:,1); %) Pull out times
    allcell=spkfiles{thismodel}(:,2); % Pull out cell populations
    allpops=spkfiles{thismodel}(:,3); % Pull out cell populations
    matches=(allpops==popinds(j)); % Find which cells belong to this population
    matches=(matches & rand(size(matches))<fractoplot); % Downsample randomly
    cellmatches=cellids(cellpops==popinds(j),1); % Cells that match
    scatter(alltime(matches)-plotoffset,allcell(matches),3,popcolors(j,:),'filled','marker','o'); hold on
    text(xmax+0.2,mean(cellmatches)+labeloffset-(TCRoffset*(any(popinds(j)==[TCR,IL6]))),popnames{j},'color',[0 0 0],'verticalalignment','middle','fontsize',8)
    if j~=length(popinds), plot([0 xmax+0.1],max(cellmatches)*[1 1],'k','clipping','off'), end
end
xlim([0 xmax])
ylabel('Cell number')
xlabel('Time (s)')
box on



%% Plot time series
ylims=40;
squish=0.4; % It's arbitrary units so I can squish it however I want
colors=vectocolor(1:nfiles+1)/1.3;
for f=1:nfiles
    cksubplot([2 8], [2 f], [90,65], [12,-25], [17,0], [0,0])
    thislfp=lfpfiles{f}(:,2);
    thislfp=thislfp-mean(thislfp);
    thislfp=thislfp/1000;
    plot(lfpfiles{f}(:,1)-plotoffset,thislfp*squish,'color',colors(f,:),'linewidth',2) % This is an offset because the beginning of the LFP looks funky
    xlim([0 xmax])
    ylim(ylims*[-1 1])
    if f==3, ylabel('                                      Layer 2/3 LFP voltage (arbitrary units)'), end
    if f~=nfiles, set(gca,'xticklabel',[]), else xlabel('Time (s)'), end
    title(fullkw{f})
end





%% Plot firing rates and Fano factors
disp('Plotting firing rates and Fano factors...')

% Define different ways to group the cells
disp('  Defining population -> group...')
groupings={popinds,[E2 E4 E5B E5R E6 TCR],[I2 I4 I5 I6 TRN],[IL2 IL4 IL5 IL6]};
ngroupings=length(groupings);
groupingnames={'All','Excitatory','Fast-spiking','Low-threshold','Layer 2/3','Layer 4','Layer 5','Layer 6','Thalamus'};
pop2group=cell(max(popinds),1);
for p=1:max(popinds)
    for g=1:ngroupings
        if any(groupings{g}==p)
            pop2group{p}=[pop2group{p} g]; % Create a way to go straight from population to group number
        end
    end
end

% Find a way of pulling out all the cell corresponding to a given group
disp('  Defining group -> cell...')
group2cell=cell(ngroupings,1);
for c=1:ncells
    thiscellid=locfiles{f}(c,1);
    thispop=locfiles{f}(c,2);
    for g=1:ngroupings
        if any(groupings{g}==thispop)
            group2cell{g}=[group2cell{g} thiscellid];
        end
    end
end
    

% Turn spike times into a raster
disp('  Binning spikes...')
maxtime=ceil(lfpfiles{1}(end,1));
binsize=min(20,maxtime); % Bin size in s
fullraster=zeros(nfiles,ncells,ceil(maxtime/binsize)+1);
for f=1:nfiles
    thesespikes=round(spkfiles{f}(:,1)/binsize)+1; % Convert from time to bin index
    thesecells=spkfiles{f}(:,2);
    nspks=length(thesespikes);
    for s=1:nspks
        fullraster(f,thesecells(s),thesespikes(s))=fullraster(f,thesecells(s),thesespikes(s))+1;
    end
end

% Select subsets of the raster to do computations on
disp('  Analyzing...')
means=zeros(nfiles,ngroupings);
varis=zeros(nfiles,ngroupings);
stds= zeros(nfiles,ngroupings);
coevs=zeros(nfiles,ngroupings);
fanos=zeros(nfiles,ngroupings);
for f=1:nfiles
    for g=1:ngroupings
        thisraster=fullraster(f,group2cell{g},:);
        means(f,g)=mean(thisraster(:))/binsize;
        varis(f,g)=var(thisraster(:))/binsize^2; % Normalize
        stds(f,g)= std(thisraster(:));
        coevs(f,g)=stds(f,g)/means(f,g);
        fanos(f,g)=varis(f,g)/means(f,g);
    end
end



% Turn spike times into a binned time series
disp('  Binning spikes properly...')
maxtime=ceil(lfpfiles{1}(end,1));
binsize=0.001; % Bin size in s -- doesn't really matter since rebinned later anyway
fullbinned=zeros(nfiles,ceil(maxtime/binsize)+1);
for f=1:nfiles
    maxspks=length(spkfiles{4}(:,1));
    nspks=length(spkfiles{f}(:,1));
    spkinds=1:nspks;
    toremove=randi(nspks,[round(1.1*(nspks-maxspks)),1]); % The 1.1 is to roughly correct for duplicates
    spkinds(toremove)=[];
    for s=spkinds
        bin=floor(spkfiles{f}(s,1)/binsize)+1;
        fullbinned(f,bin)=fullbinned(f,bin)+1;
    end
end

% Calculate the correlation coefficient for different bins
disp('  Calculating the correlation coefficient...')
maxbin=round(maxtime/2/binsize);
bins=unique(round(exp(linspace(0,log(maxbin),100))));
coeffvars=zeros(length(bins),4);
for t=1:4
    for bi=1:length(bins)
        b=bins(bi);
        newlength=floor(length(fullbinned(t,:))/b);
        thisvec=zeros(newlength,1);
        for i=1:newlength
            start=(i-1)*b+1;
            finish=i*b;
            thisvec(i)=sum(fullbinned(t,start:finish));
        end
        coeffvars(bi,t)=var(thisvec)/mean(thisvec); 
    end
end




% Plot firing rates
cksubplot([2 2], [1 2], [80,40], [0,0], [11,-6], [0,0]); hold on
width=0.2;
yaxis=cell(nfiles,1);
skip=[1 (2:4)+0.3]; % Make a gap between transmitter and layer groups
for f=1:nfiles
    yaxis{f}=skip+(f-nfiles/2)*width;
    barh(yaxis{f},means(f,:),'barwidth',width,'facecolor',colors(f,:),'linestyle','none')
end
ylim([0.5 ngroupings+1.1])
set(gca,'ytick',skip)
set(gca,'yticklabel',groupingnames)
set(gca,'ydir','reverse')
xlabel('Mean firing rate (Hz)')
legend(fullkw)

% Plot coefficients of variation
cksubplot([2 2], [2 2], [80,40], [0,0], [11,-6], [1,1]); hold on
linestyles={'--','--','-','-'};
for f=1:nfiles
    smoothcoeffs=cksmooth(coeffvars(:,f),5);
    plotcoeffs=[coeffvars(1:ceil(length(bins)/2),f); smoothcoeffs((ceil(length(bins)/2)+1):end)]; % Only smooth large bin sizes
    plot(bins*binsize,plotcoeffs,'color',colors(f,:),'linewidth',2,'linestyle',linestyles{f}); hold on
end
ylim([20 1e5])
set(gca,'ytick',10.^[2:5])
set(gca,'xscale','log')
set(gca,'yscale','log')
legend(fullkw,'location','northwest')
xlabel('Bin size (s)')
ylabel('Fano factor')
