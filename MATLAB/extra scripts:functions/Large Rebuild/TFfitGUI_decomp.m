%% Decomposition of the TCfitGUI

figure(1);
set(gcf, 'Visible', 'off', 'Toolbar', 'figure', 'Color', [.8 .8 .8], ...
    'Position', [5, 5, 900 800], 'Name', 'Tuning Curves');
axis off;

ax1 = axes('Units', 'normalized', 'Position', [.06 .65 .5 .3]);
ax2 = axes('Units', 'normalized', 'Position', [.63  .36 .35 .2]);
ax3 = axes('Units', 'normalized', 'Position', [.63  .045 .35 .2]);
ax4 = axes('Units', 'normalized', 'Position', [.06  .36 .5 .2]);
ax5 = axes('Units', 'normalized', 'Position', [.06  .045 .5 .2]);
ax = [ax1 ax2 ax3 ax4 ax5];

defaultparamsBMOC = csvread('defaultparamsBMOC.csv');

%% Initial parameter values
tcnum       = 5;
alphabm     = defaultparamsBMOC(10,tcnum);
beta1bm     = defaultparamsBMOC(11,tcnum);
delta1bm    = defaultparamsBMOC(12,tcnum);
alphaoc     = defaultparamsBMOC(13,tcnum);
beta1oc     = defaultparamsBMOC(14,tcnum);
delta1oc    = defaultparamsBMOC(15,tcnum);
c21         = defaultparamsBMOC(16,tcnum);
c12         = defaultparamsBMOC(17,tcnum);
thresh      = defaultparamsBMOC(18,tcnum);

clear defaultparamsBMOC
TCnumLim   = [1 11];


%% Make text boxes and labels for parameter input

% Panel for parameter input
hp1 = uipanel('Title','Parameters',...
    'FontSize',11,'FontName','Helvetica','FontUnits','normalized',...
    'BackgroundColor',[.8 .8 .8],'Position',[.71 .66 .28 .3]);
axes('parent', hp1,'position',[0 0 1 1], 'visible', 'off');

% Label which layer the textboxes correspond to
text(.55,.96,'BM','HorizontalAlignment','center','VerticalAlignment','bottom',...
    'FontSize',12,'FontUnits','normalized','Units','normalized')
text(.85,.96,'OC','HorizontalAlignment','center','VerticalAlignment','bottom',...
    'FontSize',12,'FontUnits','normalized','Units','normalized')

% Label and textbox for alpha_BM
text(.45,.89,'alpha:','HorizontalAlignment','right','VerticalAlignment','middle',...
    'FontSize',12,'FontUnits','normalized','Units','normalized')
handles.alphabm = uicontrol('Parent',hp1,'Style','edit', 'String', alphabm,...
    'FontSize',12,'FontUnits','normalized','Units','normalized','Position', [.45 .83 .2 .13]);
% Label and textbox for beta_BM
text(.45,.74,'beta:','HorizontalAlignment','right','VerticalAlignment','middle',...
    'FontSize',12,'FontUnits','normalized','Units','normalized')
handles.betabm = uicontrol('Parent',hp1,'Style','edit', 'String', beta1bm,...
    'FontSize',12,'Units', 'normalized','FontUnits','normalized','Position', [.45 .68 .2 .13],'enable','off');
% Label and textbox for delta_BM
text(.45,.59,'delta:','HorizontalAlignment','right','VerticalAlignment','middle',...
    'FontSize',12,'FontUnits','normalized','Units','normalized')
handles.deltabm = uicontrol('Parent',hp1,'Style','edit', 'String', delta1bm,...
    'FontSize',12,'Units', 'normalized','FontUnits','normalized','Position', [.45 .53 .2 .13],'enable','off');

% Label and textbox for alpha_OC
handles.alphaoc = uicontrol('Parent',hp1,'Style','edit', 'String', alphaoc,...
    'FontSize',12,'Units', 'normalized','FontUnits','normalized','Position', [.76 .83 .2 .13]);
% Label and textbox for beta_OC
handles.betaoc = uicontrol('Parent',hp1,'Style','edit', 'String', beta1oc,...
    'FontSize',12,'Units', 'normalized','FontUnits','normalized','Position', [.76 .68 .2 .13]);
% Label and textbox for delta_OC
handles.deltaoc = uicontrol('Parent',hp1,'Style','edit', 'String', delta1oc,...
    'FontSize',12,'Units', 'normalized','FontUnits','normalized','Position', [.76 .53 .2 .13]);

% Label and textbox for c12
text(.45,.45,'c12:','HorizontalAlignment','right','VerticalAlignment','middle',...
    'FontSize',12,'FontUnits','normalized','Units','normalized')
handles.c12 = uicontrol('Parent',hp1,'Style','edit', 'String', c12,...
    'FontSize',12,'Units', 'normalized','FontUnits','normalized','Position', [.45 .38 .2 .13]);
% Label and texbox for c21
text(.76,.45,'c21:','HorizontalAlignment','right','VerticalAlignment','middle',...
    'FontSize',12,'FontUnits','normalized','Units','normalized')
handles.c21 = uicontrol('Parent',hp1,'Style','edit', 'String', c21,...
    'FontSize',12,'Units', 'normalized','FontUnits','normalized','Position', [.76 .38 .2 .13]);

% Label and textbox for threshold
text(.05,.83,'Threshold:','HorizontalAlignment','left',...
    'VerticalAlignment','middle','Units','normalized','FontSize',12,'FontUnits','normalized')
handles.thresh = uicontrol('Parent', hp1, 'Style','edit', 'String', thresh,...
    'FontSize',12,'Units', 'normalized', 'HorizontalAlignment', 'center',...
    'Position', [.05 .67 .2 .13],'FontUnits','normalized');

% Pushbutton controls for saving as default parameters and choosing c12
% value
handles.buttonSaveAsDefault = uicontrol('Parent', hp1, 'Style', 'pushbutton', 'String', {'Save as Default'},...
    'TooltipString', ['Save the current parameters' char(10)...
    'as the default parameters in' char(10) 'defaultparamsBMOC.csv'],...
    'FontUnits','normalized','Units', 'normalized',  'Position', [.5 .2 .4 .18]);
handles.buttonChoosec12 = uicontrol('Parent', hp1, 'Style', 'pushbutton', 'String', 'Choose c12',...
    'FontUnits','normalized','Units', 'normalized',...
    'TooltipString', ['Set value of c12 to make spontaneous' char(10) 'amplitude equal to half the threshold.'],...
    'Position', [0 .37 .35 .15]);

% Tuning curve slider, textbox, and label
% axes('position',[.8 .62 .05 .02],'Visible','off');
text(.33,.13,'Tuning Curve:','HorizontalAlignment','right','VerticalAlignment','baseline',...
    'FontSize',12,'FontUnits','normalized','Units','normalized')
handles.textTC = uicontrol('Style','edit', 'String', tcnum,...
    'FontSize',12,'FontUnits','normalized','Units','normalized','Position', [.805 .69 .046 .03]);
handles.sliderTC = uicontrol('Style', 'slider', 'Min', 1, 'Max', 11, 'Value', tcnum, 'SliderStep', [.1 .4], ...
    'Min', TCnumLim(1), 'Max', TCnumLim(2), 'Units', 'normalized', 'Position', [.72 .66 .15 .025]);


% Create panel for toggle controls
hp2 = uipanel('Title','Toggles',...
    'FontSize',11,'FontName','Helvetica','FontUnits','normalized',...
    'BackgroundColor',[.8 .8 .8],'Position',[.57 .66 .13 .3]);
axes('Units','normalized','Position',get(hp2,'Position'))
axis off

% Toggle controls for middle ear filter, frequency scaling, setting organ
% of corti to be the treshold, and using a single oscillator network
handles.toggleMEF = uicontrol('Parent',hp2,'Style','togglebutton', 'String','MEF',...
    'TooltipString', 'Turn ON/OFF Middle Ear Filter',...
    'Value', 1,'FontUnits','normalized','Units', 'normalized', 'Position', [.03 0.8 0.95 0.15]);
handles.toggleFS = uicontrol('Parent',hp2,'Style','togglebutton', 'String','Freq Scaling',...
    'TooltipString', 'Turn ON/OFF Frequency Scaling',...
    'Value', 0,'FontUnits','normalized','Units', 'normalized', 'Position', [.03 .55 .95 .15]);
handles.toggleThresh = uicontrol('Parent',hp2,'Style','togglebutton', 'String','Thresh: OC',...
    'TooltipString', ['ON = use OC amplitude as threshold ' char(10) 'OFF = use BM amplitude as threshold'],...
    'Value', 1,'FontUnits','normalized','Units', 'normalized', 'Position', [.03 .3 .95 .15]);
handles.toggleSingle = uicontrol('Parent',hp2,'Style','togglebutton', 'String','Two Layer',...
    'TooltipString', ['ON = two layer model' char(10) 'OFF = single layer model'],...
    'Value', 1,'FontUnits','normalized','Units', 'normalized', 'Position', [.03 .05 .95 .15]);


%% Set callback functions
handlesTC           = [handles.sliderTC handles.textTC];
handles.textGroup   = [handles.alphabm handles.betabm handles.deltabm ...
    handles.alphaoc handles.betaoc handles.deltaoc ...
    handles.c21 handles.c12 handles.thresh];
handles.toggleGroup = [handles.toggleMEF handles.toggleFS handles.toggleThresh handles.toggleSingle];

set(handlesTC,                  'Callback', {@Callback_TC, handles, ax});
set(handles.textGroup,          'Callback', {@text_Callback,   handles, ax});
set(handles.buttonSaveAsDefault,'Callback', {@saveAsDef_Callback, handles});
set(handles.buttonChoosec12,    'Callback', {@choosec12_Callback, handles, ax});
set(handles.toggleGroup,        'Callback', {@toggle_Callback,handles,ax});

hold off

toggleValues = {1, 0, 1, 1};

% Load the Empirical Tuning Curve data.
load('JorisTCdata.mat');

% Define Reference Value
Fref    = 0.00002;  % Reference pressure in Pa.

tcnum = round(tcnum);


cf            = JorisTCdata{tcnum}(2,1);     % Center Frequency CF.
thresholdAtCF = JorisTCdata{tcnum}(1,1);
f0            = JorisTCdata{tcnum}(2,2:end); % Input/Stimulus Frequency
thresholds    = JorisTCdata{tcnum}(1,2:end);

idxcf         = f0 == cf;

f        = cf;
fsquared = f^2;

omega        =  2*pi*(f0 - f);
omegasquared = (2*pi*(f0 - f)).^2;
omega2       = omegasquared/fsquared;

MEFdBgains   = MEFdBgain(f0);  % In dB

pAt0dB = dB2Pa( MEFdBgain( [f-1 f f+1] ));
pAt0dB = pAt0dB(2);


rbm = thresh/c21*sqrt((alphaoc + beta1oc*thresh^2)^2 + (omega + delta1oc*thresh^2).^2);
sinPsioc = thresh*(omega + delta1oc*thresh^2)/(c21*rbm);
cosPsioc = sqrt(1-sinPsioc.^2);
F = sqrt((alphabm*rbm + c12*thresh*cosPsioc).^2 + (omega.*rbm + c12*thresh*sinPsioc).^2);    % Ji Chul's calculation for bidirectional coupling

FStable   = NaN(1,length(omega));
FUnstable = NaN(3,length(omega));

for j=1:length(omega)
    for i=1:size(F,1)
        
        F = sqrt((alphabm*rbm + c12*thresh*cosPsioc).^2 + (omega.*rbm + c12*thresh*sinPsioc).^2);
        % not frequency scaling
        %%%%%%%
        %         [~,r2,~,~,stable] = rStarCochlea(alphabm,alphaoc,beta1oc,delta1oc,F(i,j),c21,c12,omega(j),1);
        alphaBM = alphabm;
        alphaOC = alphaoc;
        beta1OC = beta1oc;
        delta1OC = delta1oc;
        F = F(i,j);
        Omega = omega(j);
        All = 1;
        
        abm = alphaBM; aoc = alphaOC; b1oc = beta1OC; d1oc = delta1OC;
        W = Omega;
        
        % Get rOC*
        roc = sqrt(roots([ ...
            power(power(b1oc,2) + power(d1oc,2),2)*(power(abm,2) + power(W,2)), 2*(power(b1oc,2) + power(d1oc,2))*(-(abm*b1oc*c12*c21) + ...
            2*power(abm,2)*(aoc*b1oc + d1oc*W) + W*(c12*c21*d1oc + ...
            2*W*(aoc*b1oc + d1oc*W))), -2*abm*c12*c21*(aoc*(3*power(b1oc,2) + power(d1oc,2)) + ...
            2*b1oc*d1oc*W) + 4*aoc*b1oc*d1oc*W*(c12*c21 + 2*power(W,2)) + ...
            2*power(abm,2)*(power(aoc,2)*(3*power(b1oc,2) + power(d1oc,2)) + ...
            4*aoc*b1oc*d1oc*W + (power(b1oc,2) + 3*power(d1oc,2))*power(W,2)) + ...
            power(b1oc,2)*(power(c12,2)*power(c21,2) + 6*power(aoc,2)*power(W,2) + 2*c12*c21*power(W,2) + 2*power(W,4)) + ...
            power(d1oc,2)*(power(c12,2)*power(c21,2) + 2*power(aoc,2)*power(W,2) + 6*c12*c21*power(W,2) + 6*power(W,4)),...
            -(power(b1oc,2)*power(c21,2)*power(F,2)) - power(c21,2)*power(d1oc,2)*power(F,2) + ...
            2*power(c12,2)*power(c21,2)*d1oc*W + 4*power(aoc,3)*b1oc*power(W,2) + ...
            6*c12*c21*d1oc*power(W,3) + 4*d1oc*power(W,5) + 4*power(abm,2)*(aoc*b1oc + d1oc*W)*(power(aoc,2) + power(W,2)) + ...
            2*power(aoc,2)*d1oc*W*(c12*c21 + 2*power(W,2)) - ...
            2*abm*c12*c21*(3*power(aoc,2)*b1oc + 2*aoc*d1oc*W + ...
            b1oc*power(W,2)) + 2*aoc*b1oc*(power(c12,2)*power(c21,2) + ...
            2*c12*c21*power(W,2) + 2*power(W,4)), -2*aoc*b1oc*power(c21,2)*power(F,2) + power(aoc,4)*power(W,2) - ...
            2*abm*aoc*c12*c21*(power(aoc,2) + power(W,2)) + power(abm,2)*power(power(aoc,2) + power(W,2),2) + ...
            power(aoc,2)*(power(c12,2)*power(c21,2) + 2*c12*c21*power(W,2) + ...
            2*power(W,4)) + W*(2*c12*c21*power(W,3) + power(W,5) + ...
            power(c21,2)*(-2*d1oc*power(F,2) + power(c12,2)*W)), -(power(c21,2)*power(F,2)*(power(aoc,2) + power(W,2)))]));
        
        roc = roc(abs(imag(roc)) < eps('single')); % take only real roots
        roc = roc(abs(roc) > 0); % take only positive roots
        roc = real(roc);
        roc = sort(unique(roc),'descend'); % remove multiple roots
        
        % Get rBM*
        rbm = sqrt((aoc*roc+b1oc*roc.^3).^2 + (W*roc + d1oc*roc.^3).^2)/c21;
        ind = find(abs(imag(rbm)) < eps('single'));
        rbm = real(rbm(ind));
        roc = roc(ind);
        
        % Get psiOC*
        signPsiOC = ((W + d1oc*roc.^2)*c21 >= 0)*2 - 1; % sign of psi*oc
        psioc = signPsiOC.*acos(-(aoc*roc+b1oc*roc.^3)./(c21*rbm));
        ind = find(abs(imag(psioc)) < 10^(-5));
        psioc = real(psioc(ind));
        rbm = rbm(ind);
        roc = roc(ind);
        
        % Get psiBM*
        signPsiBM = (W + c12*roc.*sin(psioc)./rbm >= 0)*2 - 1; % sign of psi*oc
        psibm = signPsiBM.*acos(-(abm*rbm+c12*roc.*cos(psioc))/F);
        ind = find(abs(imag(psibm)) < 10^(-5));
        psibm = real(psibm(ind));
        rbm = rbm(ind);
        roc = roc(ind);
        psioc = psioc(ind);
        
        % Calculate stability type using Jacobian matrix
        stabType = zeros(size(roc));
        for n = 1:length(roc)
            J = NaN(4);
J(1,1) = abm;
J(1,2) = -F*sin(psibm(n));
J(1,3) = c12*cos(psioc(n));
J(1,4) = -c12*roc(n)*sin(psioc(n));
J(2,1) = (F*sin(psibm(n))-c12*roc(n)*sin(psioc(n)))/(rbm(n))^2;
J(2,2) = -F*cos(psibm(n))/rbm(n);
J(2,3) = c12*sin(psioc(n))/rbm(n);
J(2,4) = c12*roc(n)*cos(psioc(n))/rbm(n);
J(3,1) = c21*cos(psioc(n));
J(3,2) = 0;
J(3,3) = aoc + 3*b1oc*roc(n)^2;
J(3,4) = -c21*rbm(n)*sin(psioc(n));
% J(4,1) = -c21*sin(psioc)/roc;
% J(4,2) = 0;
% J(4,3) = 2*d1oc*roc + c21*rbm*sin(psioc)/roc^2;
% J(4,4) = -c21*rbm*cos(psioc)/roc;
J(4,1) = -c21*sin(psioc(n))/roc(n) - (F*sin(psibm(n))-c12*roc(n)*sin(psioc(n)))/rbm(n)^2;
J(4,2) = F*cos(psibm(n))/rbm(n);
J(4,3) = 2*d1oc*roc(n) + (c21*rbm(n)/roc(n)^2 - c12/rbm(n))*sin(psioc(n));
J(4,4) = -(c21*rbm(n)/roc(n) + c12*roc(n)/rbm(n))*cos(psioc(n));
            
            ev = eig(J);
            if isreal(ev) && all(ev < 0)
                stabType(n) = 4; % stable node
            elseif all(real(ev) < 0)
                stabType(n) = 3; % stable spiral
            elseif isreal(ev) && all(ev > 0)
                stabType(n) = 2; % unstable node
            elseif all(real(ev) > 0)
                stabType(n) = 1; % unstable spiral
            end % saddle pt otherwise
        end
        stable = (stabType >= 3); % 1 = stable, 0 = unstable
        
        % Prepare output
        if All % both stable and unstable solutions
          rStarBM = rbm;
          r2 = roc;
          psiStarBM = psibm;
          psiStarOC = psioc;
        else % only stable solutions
          indStab = find(stable);
          rStarBM = rbm(indStab);
          r2 = roc(indStab);
          psiStarBM = psibm(indStab);
          psiStarOC = psioc(indStab);
          stable = stable(indStab);
          stabType = stabType(indStab);
        end
        
        %%%%%%%
        [~,index] = min(abs(thresh-r2));
        %%%%%%%
        FStable(j)     = F;
        %%%%%%%
    end
end

% middle ear filtering

LstimStable   = 20*log10(FStable/Fref) - MEFdBgains;
LstimUnstable = 20*log10(FUnstable/Fref) - repmat(MEFdBgains,3,1);


semilogx(ax(1), f0, thresholds, 'b.-', 'Linewidth', 2);  hold(ax(1), 'on');
semilogx(ax(1), f0, LstimStable,   'Color', [.2 .8 .2], 'LineWidth', 2, 'Marker', '.');
semilogx(ax(1), f0, LstimUnstable, 'Color', [.8 .2 .8], 'LineWidth', 2, 'Marker', '.');

meanSquareError = mean((thresholds-LstimStable).^2);

xlim(ax(1), [50 30000]);
ylim(ax(1), [0 90]);
grid(ax(1), 'on');
ylabel(ax(1),'Threshold (dB SPL)','FontUnits','normalized');
xlabel(ax(1),'Forcing frequency (Hz)','FontUnits','normalized');
title(ax(1), ['TC ' num2str(tcnum) ' Fit  with  Tip = (' num2str(f) ', ' num2str(LstimStable(idxcf)) '), RMSE = ' num2str(sqrt(meanSquareError)) ' dB'],...
    'FontUnits','normalized');
set(ax(1), 'FontUnits','normalized');
hold(ax(1), 'off');


%% Amplitude Curves Method 2
% Compute rocs by specifying Fs then use then use the rocs to compute
% rbms

F = logspace(log10(Fref),log10(20),201);

rBM = [];
rOC = [];
F_rBM_at_rOCthreshold = [];

clear f0;
n = 0;

threef0 = [f f/2 2*f];

for f0 = threef0
    n = n+1;
    
    rstarbm = NaN(size(F));
    rstaroc = NaN(size(F));
    for nF = 1:length(F)
        % coupled oscillator model
        % not frequency scaling
        [rbmn, rocn] = rStarCochlea(alphabm, alphaoc, beta1oc, delta1oc, ...
            F(nF), c21, c12, 2*pi*(f-f0));
        rstarbm(nF) = rbmn(1);
        rstaroc(nF) = rocn(1);
    end
    rBM = [rBM; rstarbm];
    rOC = [rOC; rstaroc];
end

clear rstaroc rstarbm

%% Plot Amplitude Curves

F = Pa2dB(F); % For a log-log plot.
%F_rBM_at_rOCthreshold(:,1) = Pa2dB(F_rBM_at_rOCthreshold(:,1));

% coupled oscillator model
% oc==thresh
% afferent and efferent connections are both zero
semilogy(ax(2), F, rOC(1,:), 'g', F, rOC(2,:), 'b', F, rOC(3,:), 'r', [Fref 120], [thresh thresh], 'k--', 'Linewidth', 2); %, [0 max(F)] , [rbm rbm], 'k--');
semilogy(ax(3), F, rBM(1,:), 'g', F, rBM(2,:), 'b', F, rBM(3,:), 'r', 'Linewidth', 2);

xlabel(ax(2), 'Forcing F (dB SPL)','FontUnits','normalized');
ylabel(ax(2), 'r_{OC}^*','FontUnits','normalized');
title(ax(2), 'Steady State Amplitude Curve for OC Oscillator.',...
    'FontUnits','normalized');
set(ax(2), 'XGrid', 'on', 'YGrid', 'on', 'YMinorGrid', 'off');
set(ax(2), 'FontUnits','normalized');
ax2ylim = get(ax(2), 'YLim');
ylim(ax(2), [rOC(1,1)/10 ax2ylim(2)]);


grid(ax(3), 'on');

xlabel(ax(3), 'Forcing F (dB SPL)','FontUnits','normalized');
ylabel(ax(3), 'r_{BM}^*','FontUnits','normalized');
title(ax(3), 'Steady State Amplitude Curve for BM Oscillator.',...
    'FontUnits','normalized');
set(ax(3), 'XGrid', 'on', 'YGrid', 'on', 'YMinorGrid', 'off');
set(ax(3), 'FontUnits','normalized');
ax3ylim = get(ax(3), 'YLim');
ylim(ax(3), [rBM(1,1)/10 ax3ylim(2)]);

%% OC & BM Compression Curves

f0sOrig = logspace(log10(50),log10(30000),201)';
temp    = f0sOrig - f;
[~,ind] = min(abs(temp));
f0s     = f0sOrig - temp(ind);
F       = dB2Pa([120 100 80 60 40 20 0]);
W = 2*pi*(f-f0s);
BMcomp = zeros(length(W),length(F));
OCcomp = zeros(length(W),length(F));

for nF = 1:length(F)
    for nW = 1:length(W)
        % coupled oscillator model
        [rbmn, rocn] = rStarCochlea(alphabm, alphaoc, beta1oc, delta1oc, ...
            F(nF), c21, c12, W(nW));
            if isempty(rbmn)                
                BMcomp(nW,nF) = NaN;
                OCcomp(nW,nF) = NaN;
            else
                BMcomp(nW,nF) = rbmn(1);
                OCcomp(nW,nF) = rocn(1);
            end
    end
end

% using OC as threshold
loglog(ax(5), f0s(f0s>0), BMcomp(f0s>0,:), 'Linewidth', 2);
loglog(ax(4), f0s(f0s>0), OCcomp(f0s>0,:), 'Linewidth', 2); hold(ax(4), 'on');
loglog(ax(4), [min(f0sOrig) max(f0sOrig)], [thresh thresh], 'k--', 'Linewidth', 2); hold(ax(4), 'off');

xlim(ax(5), [min(f0sOrig) max(f0sOrig)]);
set(ax(5), 'XGrid', 'on', 'YGrid', 'on', 'YMinorGrid', 'off');
xlabel(ax(5), 'Forcing frequency (Hz)','FontUnits','normalized');
ylabel(ax(5), 'r_{BM}^*','FontUnits','normalized');
title(ax(5), 'BM Compression Curves','FontUnits','normalized');
set(ax(5), 'FontUnits','normalized');

xlim(ax(4), [min(f0sOrig) max(f0sOrig)]);

ylim(ax(4), [thresh/10 2*max(rOC(1,:))] );
set(ax(4), 'XGrid', 'on', 'YGrid', 'on', 'YMinorGrid', 'off');
xlabel(ax(4), 'Forcing frequency (Hz)','FontUnits','normalized');
ylabel(ax(4), 'r_{OC}^*','FontUnits','normalized');
title(ax(4), 'OC Compression Curves','FontUnits','normalized');
set(ax(4), 'FontUnits','normalized');

%% Decide what needs to be plotted based on toggle button values

stringL1 = {'120','100','80','60','40','20','0dB','thresh'};
stringL2 = {'f0 = CF = f','f0 = f/2','f0 = 2f'};

hl2 = legend(ax(2),stringL2,'Orientation', 'horizontal','fontsize', 12,...
    'units', 'normalized');
hl4 = legend(ax(4),stringL1,'Orientation', 'horizontal','fontsize', 12,...
    'units', 'normalized');

% Set position, orientation of legends
temp2 = get(hl2, 'Position'); hl2size = temp2(3:4);
temp4 = get(hl4, 'Position'); hl4size = temp4(3:4);

set(hl4, 'Position', [.02 .28 hl4size(:)'])
set(hl2, 'Position', [.7 .28 hl2size(:)'])


movegui(gcf,'center');
set(gcf, 'Visible', 'on');
