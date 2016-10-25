%% FIG7GRANGER
% This code plots Fig. 7 for the paper: Granger causality between
% different layers of the cortex.
% Version: 2012dec29 by Cliff Kerr (cliffk@neurosim.downstate.edu)

%% Load data
kw={'wh','tc','bg','pa'};
fullkw={'WN','TC','BG','PD'};

nfiles=length(kw);
disp('Loading data...')
base='../output/output'; % Base file name -- needs to match runsim
lfpfiles=cell(nfiles,1);
for f=1:nfiles
    fprintf('  File %i of %i...\n',f,nfiles)
    lfpfiles{f}=dlmread([base '-' kw{f} '-lfp.txt']);
    lfpfiles{f}(:,1)=lfpfiles{f}(:,1)/1000; % Convert from ms to s
end

%% Trim data
disp('Trimming data...')
trimends=[3 3]; % How many seconds of data to trim off the ends
for f=1:nfiles
    maxtime=lfpfiles{f}(end,1);
    start=trimends(1); 
    finish=maxtime-trimends(2);
    lfplimits=(lfpfiles{f}(:,1)>start & lfpfiles{f}(:,1)<finish);
    lfpfiles{f}=lfpfiles{f}(lfplimits,:);
    lfpfiles{f}(:,1)=lfpfiles{f}(:,1)-trimends(1);
end


%% Calculate all causalities
disp('Calculating causality...')
nlayers=4; % Number of cortical layers: 2/3, 4, 5, 6
fs=200; % Sampling rate in Hz
order=12; % Polynomial order for autoregressive fit
maxfreq=50; % Maximum frequency in Hz to plot
allcausality=cell(nlayers,nlayers);

for inlayer=1:nlayers
    for outlayer=1:nlayers
        if inlayer~=outlayer
            fprintf ('  %i/%i of %i/%i\n',inlayer,outlayer,nlayers,nlayers)

            causality=cell(nfiles,5);
            for f=1:nfiles
                thisdata1=lfpfiles{f}(:,inlayer+1); % +1 since first is time
                thisdata2=lfpfiles{f}(:,outlayer+1);
                thisdata1=thisdata1-mean(thisdata1);
                thisdata2=thisdata2-mean(thisdata2);
                thisdata=[thisdata1';thisdata2'];
                [causality{f,:}]=pwcausalr(thisdata,1,length(thisdata),order,fs,maxfreq); % Calculate causality
            end
            allcausality{inlayer,outlayer}=causality;
        end
    end
end




%% Create figure out of selected causalities
disp('Plotting...')
h=figure('position',[200 200 700 500]);

% Plot labels
ax=axes('position',[0 0 1 1],'units','normalized'); hold on
text(0.02,0.97,'A','fontsize',16,'fontweight','bold')
text(0.52,0.97,'B','fontsize',16,'fontweight','bold')
text(0.02,0.49,'C','fontsize',16,'fontweight','bold')
text(0.52,0.49,'D','fontsize',16,'fontweight','bold')
xlim([0 1])
ylim([0 1])
set(ax,'visible','off')

% Plot content
pairs={[2 1],[1 3],[2 3],[4 1]};
labels={'2/3','4','5','6'};
linestyles={'--','--','-','-'};
colors=vectocolor(1:nfiles+1)/1.3;
xpos=[1 2 1 2];
ypos=[1 1 2 2];
for p=1:length(pairs)
    cksubplot([2 2],[xpos(p) ypos(p)],[95 75],[5 5],[7 -3],[0 0])
    for f=1:nfiles
        plot(0:maxfreq,allcausality{pairs{p}(1),pairs{p}(2)}{f,3},'color',colors(f,:),'linewidth',2,'linestyle',linestyles{f}); hold on
        ylim([0 0.8])
    end
    if xpos(p)==1, ylabel('Granger causality')
    else set(gca,'yticklabel',[])
    end
    if ypos(p)==1, set(gca,'xticklabel',[])
    else xlabel('Frequency (Hz)')
    end
    title(['Layer ' labels{pairs{p}(1)} ' \rightarrow ' labels{pairs{p}(2)}])
    if p==1, legend(fullkw,'location','northwest'), end
end
