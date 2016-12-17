function generate_fig(Features_opt, plot_flag)

% function plotting(Features_opt, data_type)
% This contains functions for plotting output from model
% INPUTS
% Features_opt - A vector of features (firing rates and frequency)
%                obtained from the model.
% plot_flag    - If plot_flag = 1 then it plots figures from 
%                full_long_weight20_pooled.mat. If plot_flag = 2 then it
%                plots figures from zerowsc_long_weight20_pooled.mat. 

% Set plot layout parameters
fontsize = 5;
labelsize = 6;
annotation_horz_pos = -1.6;
annotation_vert_pos = 0.5;
annotate_fontsize = 8;

% extract variables
minSTN  = cell2mat(Features_opt(1));
meanSTN = cell2mat(Features_opt(2));
maxSTN  = cell2mat(Features_opt(3));
minGPe  = cell2mat(Features_opt(4));
meanGPe = cell2mat(Features_opt(5));
maxGPe  = cell2mat(Features_opt(6)); 
freq    = cell2mat(Features_opt(7)); 
x1      = cell2mat(Features_opt(8));  %STN firing rate of intact model
y1      = cell2mat(Features_opt(9));  %GPe firing rate of intact model
x2      = cell2mat(Features_opt(10)); %STN firing rate of model with wgs =0
y2      = cell2mat(Features_opt(11)); %GPe firing rate of modelwith wgs = 0
x3      = cell2mat(Features_opt(12)); %STN firing rate of model with wsg = 0
y3      = cell2mat(Features_opt(13)); %GPe firing rate of model with wsg = 0
x4      = cell2mat(Features_opt(14)); %STN firing rate of model with wcs = 0
y4      = cell2mat(Features_opt(15)); %GPe firing rate of model with wcs = 0
x5      = cell2mat(Features_opt(18)); %str = 0
y5      = cell2mat(Features_opt(19)); %str = 0
x6      = cell2mat(Features_opt(16)); %STN firing rate of model with wsc = 0
y6      = cell2mat(Features_opt(17)); %GPe firing rate of model with wsc = 0


%% 
    % firing rates of single neurons
    experimental = [5, 65, 125, 45, 100, 155, 14];  
    simulation = [minSTN, meanSTN, maxSTN, minGPe, meanGPe, maxGPe, freq];
    
if plot_flag == 1;    
    hFig = figure;
    set(hFig, 'Position', [0 0 700 1000])   
    plot_flag_vec = [1,3,5,7,9,11,13];
elseif plot_flag == 2; 
    plot_flag_vec = [2,4,6,8,10,12,14];
    hold on;
end

subplot(7,2, plot_flag_vec(7))
plot(experimental, 'x',...
    'LineStyle','none',...
    'LineWidth', 2, ...
    'MarkerSize',6,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[0, 0, 1])

set(gca, 'Units', 'Normalized');
Labels = {'minSTN','meanSTN','maxSTN','minGPe','meanGPe','maxGPe','freq'};
xticks = 1:7; 
set(gca, ...
    'Box'         , 'on'                        , ...
    'LooseInset'  , get(gca, 'TightInset') * 1.5 , ...
    'TickDir'     , 'in'                         , ...
    'TickLength'  , [.02 .02]                    , ...
    'LineWidth'   , 1                            , ...
    'XMinorTick'  , 'off'                        , ...
    'YMinorTick'  , 'off'                        , ...
    'XTick'       , xticks                       , ...
    'XTickLabel'  , Labels                       , ...
    'FontSize'    , fontsize                         );
rotateXLabels(gca(), 90)

hold on;
subplot(7,2,plot_flag_vec(7))
plot(simulation,'o',...
    'LineStyle','none',...
    'MarkerSize',5,...
    'MarkerEdgeColor','b',...
    'MarkerFaceColor',[1, 0, 0])

grid on
if plot_flag == 1;

h_legend = legend('experimental','simulation');
lh=findall(gcf,'tag','legend');
lp=get(lh,'position');
set(lh,'position',[0.40,0.185,lp(3:4)]);

set(h_legend,'FontSize',fontsize);
ylabel({'Firing rate (spk/s) or' ;'oscillation frequency (Hz)'},'fontsize',fontsize)
hold off
end

%% 

subplot(7,2,plot_flag_vec(1))
plot(x1,y1)
if plot_flag == 1;
ylabel('Firing rate (spk/s)','fontsize',fontsize)
end
if plot_flag == 1;
h_legend = legend('STN','GP', 'E', 'I');
set(h_legend,'FontSize',fontsize);
end

% Headers for set of subplots
if plot_flag == 1;
    descr = {'Resonance Model'};
    ax = gca;
    axes(ax) 
    text(0.3,250,descr,'fontsize',10);
elseif plot_flag == 2;
    descr = {'Feedback Model'};
    ax = gca;
    axes(ax) 
    text(0.3,370,descr,'fontsize',10);
end

set(gca, 'Units', 'Normalized');
set(gca, ...
    'Box'         , 'on'                        , ...
    'LooseInset'  , get(gca, 'TightInset') * 1.5 , ...
    'TickDir'     , 'in'                         , ...
    'TickLength'  , [.02 .02]                    , ...
    'LineWidth'   , 1                            , ...
    'FontSize'    , labelsize                           );

descr = {'Intact'};
ax = gca;
axes(ax) 
h = text(annotation_horz_pos,annotation_vert_pos,descr,'fontsize',annotate_fontsize);
set(h, 'rotation', 90)

%% 
subplot(7,2,plot_flag_vec(2))
plot(x2,y2)
if plot_flag == 1;
ylabel('Firing rate (spk/s)','fontsize',fontsize)
end

set(gca, 'Units', 'Normalized');
set(gca, ...
    'Box'         , 'on'                        , ...
    'LooseInset'  , get(gca, 'TightInset') * 1.5 , ...
    'TickDir'     , 'in'                         , ...
    'TickLength'  , [.02 .02]                    , ...
    'LineWidth'   , 1                            , ...
    'FontSize'    , labelsize                               );

descr = {'W_{GS}=0'};
ax = gca;
axes(ax) 
h = text(annotation_horz_pos,annotation_vert_pos,0,descr,'fontsize',annotate_fontsize);
set(h, 'rotation', 90)

subplot(7,2,plot_flag_vec(3))
plot(x3,y3)
if plot_flag == 1;
ylabel('Firing rate (spk/s)','fontsize',fontsize)
end

set(gca, 'Units', 'Normalized');
set(gca, ...
    'Box'         , 'on'                        , ...
    'LooseInset'  , get(gca, 'TightInset') * 1.5 , ...
    'TickDir'     , 'in'                         , ...
    'TickLength'  , [.02 .02]                    , ...
    'LineWidth'   , 1                            , ...
    'FontSize'    , labelsize                               );

descr = {'W_{SG}=0'};
ax = gca;
axes(ax) 
h = text(annotation_horz_pos,annotation_vert_pos,0,descr,'fontsize',annotate_fontsize);
set(h, 'rotation', 90)

%%
subplot(7,2,plot_flag_vec(4))
plot(x4,y4)
if plot_flag == 1;
ylabel('Firing rate (spk/s)','fontsize',fontsize)
end

set(gca, 'Units', 'Normalized');
set(gca, ...
    'Box'         , 'on'                        , ...
    'LooseInset'  , get(gca, 'TightInset') * 1.5 , ...
    'TickDir'     , 'in'                         , ...
    'TickLength'  , [.02 .02]                    , ...
    'LineWidth'   , 1                            , ...
    'FontSize'    , labelsize                               );

descr = {'W_{CS}=0'};
ax = gca;
axes(ax) 
h = text(annotation_horz_pos,annotation_vert_pos,0,descr,'fontsize',annotate_fontsize);
set(h, 'rotation', 90)

%% 
subplot(7,2,plot_flag_vec(5))
plot(x5,y5)
if plot_flag == 1;
ylabel('Firing rate (spk/s)','fontsize',fontsize)
end

set(gca, 'Units', 'Normalized');
set(gca, ...
    'Box'         , 'on'                        , ...
    'LooseInset'  , get(gca, 'TightInset') * 1.5 , ...
    'TickDir'     , 'in'                         , ...
    'TickLength'  , [.02 .02]                    , ...
    'LineWidth'   , 1                            , ...
    'FontSize'    , labelsize                               );

descr = {'Str = 0'};
ax = gca;
axes(ax) 
h = text(annotation_horz_pos,annotation_vert_pos,0,descr,'fontsize',annotate_fontsize);
set(h, 'rotation', 90)

%% 
subplot(7,2,plot_flag_vec(6))
plot(x6,y6)
if plot_flag == 1;
ylabel('Firing rate (spk/s)','fontsize',fontsize)
end
xlabel('Time (s)','fontsize',fontsize)
xlabh = gca();

if plot_flag == 1;
    ylimits = ylim; 
set(get(xlabh, 'XLabel'), 'Position', [0.5 -ylimits(2)/6 0])
elseif plot_flag == 2;
    ylimits = ylim; 
set(get(xlabh, 'XLabel'), 'Position', [0.5 -ylimits(2)/6 0])
end

set(gca, 'Units', 'Normalized');
set(gca, ...
    'Box'         , 'on'                        , ...
    'LooseInset'  , get(gca, 'TightInset') * 1.5 , ...
    'TickDir'     , 'in'                         , ...
    'TickLength'  , [.02 .02]                    , ...
    'LineWidth'   , 1                            , ...
    'FontSize'    , labelsize                               );

descr = {'w_{SC} = 0'};
ax = gca;
axes(ax) 
h = text(annotation_horz_pos,annotation_vert_pos,0,descr,'fontsize',annotate_fontsize);
set(h, 'rotation', 90)

%% 
save_handle = subplot(7,2,plot_flag_vec(6));
saveas(save_handle, 'output_figure')

end 