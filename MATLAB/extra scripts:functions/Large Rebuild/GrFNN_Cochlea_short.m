%% Decomposition of the Cochlear Model by Ed Large at al. (unpublished)
% this is the SHORT unidirectional version of the model

clear all
clc
close all

%% Utilities

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sets Matlab's current
% directory
%
% Adds all of the parent 
% folder's subdirectories
%%%%%%%%%%%%%%%%%%%%%%%%%%%
tmp = matlab.desktop.editor.getActive;
file = tmp.Filename;
clear tmp
cd(fileparts(file));
file(find(file=='/',1,'last'):end) = [];
addpath(genpath(file))
clear file
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initialize a cochlea network with 2 layers =============================

% define identity and class of the 1st layer
n1.id = 1;
n1.class = 'network';
n1.nClass = 2;

% initialize parameters for connections
n1.con = {}; % cell array containing connections targeting this network
n1.learnList = []; % indices for learned connections (used in integrator)

% parameters to define eigenvalues in the layer
n1.nFspac = 2;
n1.tick  = [];
n1.fspac = 'log';
lf   = 30;            % min freq
hf   = 10000;            % max freq
n1.N  = 831;            % number of frequency steps

% initialize the frequencies for the oscillators in the network
n1.f  = logspace(log10(lf),log10(hf),n1.N)';

% parameters for the graphical display of the network over time
n1.dStep = 50;
n1.sStep = 1;

% parameters determining the dynamics of oscillators in the network
z0 = 0;
alpha   = -412;
beta1   = 0;
beta2   = 0;
delta1  = 0;
delta2  = 0;
epsilon = 0;
n1.a  = alpha + 1i*2*pi.*n1.f;
n1.b1 = beta1 + 1i*delta1;
n1.b2 = beta2 + 1i*delta2;
n1.e = epsilon;

% initialize the state of the oscillators
n1.z0 = repmat(z0,n1.N,1);
n1.z = n1.z0;

% define identity and class of the layer
n2.id = 2;
n2.class = 'network';
n2.nClass = 2;

% initialize parameters for connections
n2.con = {}; % cell array containing connections targeting this network
n2.learnList = []; % indices for learned connections (used in integrator)

% parameters to define eigenvalues in the layer
n2.nFspac = 2;
n2.tick  = [];
n2.fspac = 'log';
lf   = 30;            % min freq
hf   = 10000;            % max freq
n2.N  = 831;            % number of frequency steps

% initialize the frequencies for the oscillators in the network
n2.f  = logspace(log10(lf),log10(hf),n2.N)';

% parameters for the graphical display of the network over time
n2.dStep = 50;
n2.sStep = 1;

% parameters determining the dynamics of oscillators in the network
z0 = 0;
alpha   = 0;
beta1   = -40816;
beta2   = 0;
delta1  = 0;
delta2  = 0;
epsilon = 0;
n2.a  = alpha + 1i*2*pi.*n2.f;
n2.b1 = beta1 + 1i*delta1;
n2.b2 = beta2 + 1i*delta2;
n2.e = epsilon;

% initialize the state of the oscillators
n2.z0 = repmat(z0,n2.N,1);
n2.z = n2.z0;

%% Make a stimulus ============================================================

% identity and class of the stimulus
s.id = 1;
s.class = 'stimulus';
s.nClass = 1; % numerical class

% parameters to generate a sinusoidal stimulus
F = dB2Pa(10); % amplitude of sinusoid in pascals
freqs = 1850; % frequency of oscillator in Hz
[NF, FREQ] = meshgrid(n1.f, freqs); % one stim freq per oscillator freq
[~, index] = min(abs(log2(NF ./ FREQ)), [], 2); % index of frequency in oscillators closer to stim freq
freqs = n1.f(index); % redefine freqs to be the closer frequency in the oscillators

% parameters for the graphical display of the stimulus over time
s.dStep = 200;
s.dispChan = 1;
s.f = [];
s.fspac = '';
s.nFspac = 0;
s.tick = [];

% DSP parameterts for the stimulus
s.ts = [0 .10];         % time spans
s.fs = 100000;         % sampling frequency
s.dt = 1/s.fs;
s.fc = num2cell(freqs);         % carrier frequencies
s.ac = num2cell(F);         % carrier amplitudes
s.sc = 0.010;       % ramp

% generate time vector of stimulus
s.t  = min(min(s.ts)):s.dt:max(max(s.ts));

% generate stim: complex sinusoid with carrier frequency and amplitude
xb = exp(1i*(2*pi.*s.fc{1,1}.*s.t));
temp = s.ac{1,1}.*xb;

% apply a ramp to the stimulus
cutoff = s.fs*s.sc;
temp = temp.';
len = size(temp,1);
scfcn = ones(len,1);
for i = 1:cutoff;
    scfcn(i) = ((i-1)/cutoff);
    scfcn(len+1-i) = ((i-1)/cutoff);
end
temp = temp.* scfcn;
s.x = temp.';

s.N = size(s.x, 1);     % number of stimulus channels

% apply middel ear filter to stimulus
tdres = 1/s.fs;
fp    = 1e3;  % prewarping frequency 1 kHz
C     = 2*pi*fp/tan(((2*pi)/2)*fp*tdres);
m11 = 1/(power(C,2) + 5.9761e+003*C + 2.5255e+007);
m12 = (-2*power(C,2) + 2*2.5255e+007);
m13 = (power(C,2) - 5.9761e+003*C + 2.5255e+007);
m14 = (power(C,2) + 5.6665e+003*C);
m15 = -2*power(C,2);
m16 = (power(C,2) - 5.6665e+003*C);
m21 = 1/(power(C,2) + 6.4255e+003*C + 1.3975e+008);
m22 = (-2*power(C,2) + 2*1.3975e+008);
m23 = (power(C,2) - 6.4255e+003*C + 1.3975e+008);
m24 = (power(C,2) + 5.8934e+003*C + 1.7926e+008);
m25 = (-2*power(C,2) + 2*1.7926e+008);
m26 = (power(C,2)-5.8934e+003*C+1.7926e+008);
m31 = 1/(power(C,2) + 2.4891e+004*C+1.2700e+009);
m32 = (-2*power(C,2) + 2*1.2700e+009);
m33 = (power(C,2) - 2.4891e+004*C + 1.2700e+009);
m34 = (3.1137e+003*C + 6.9768e+008);
m35 = 2*6.9768e+008;
m36 = (-3.1137e+003*C + 6.9768e+008);
sos = [m11*m14 m11*m15 m11*m16 1 m11*m12 m11*m13;...
    m21*m24 m21*m25 m21*m26 1 m21*m22 m21*m23;...
    m31*m34 m31*m35 m31*m36 1 m31*m32 m31*m33];
[num, den] = sos2tf(sos, 1); % gain of 1 used
s.x = filter(num, den, s.x);

%% Connect stimulus to BM =================================================


% identity and class of the connection
con.id = 1;
con.class = 'connection';
con.nClass = 3;

% obtain attributes from the source
con.source = s.id;
con.sourceClass = s.class;
con.nSourceClass = s.nClass;   % numerical class (1: stimulus, 2: network)
con.sourceAxis = s.f;
con.sourceAxisScale = s.fspac;
con.nSourceAxisScale = s.nFspac;
con.sourceAxisTick = s.tick;
con.sourceN = s.N;

% obtain attributes from the target
con.target = n1.id;
con.targetClass = n1.class;
con.nTargetClass = n1.nClass;
con.targetAxis = n1.f;
con.targetAxisScale = n1.fspac;
con.nTargetAxisScale = n1.nFspac;
con.targetAxisTick = n1.tick;
con.targetN = n1.N;

% parameters for the graphical display of the connection over time
con.type = '1freq';
con.dStep = 0;
con.sStep = 0;
con.no11 = 0;
con.tol  = .01;
con.phaseDisp = 0;
con.nType = 1; % numerical connection type
con.scale = 0;

% Connection parameters
con.F = ones(size(n1.f));
con.w = 1;
con.lambda = 0;
con.mu1 = 0;
con.mu2 = 0;
con.kappa = 0;
con.e = 0;
con.learn = 0;
con.C0 = ones(size(n1.f));
con.C = con.C0;
n1.con{con.id} = con;

%% Add connections from bm to oc ==========================================

% identity and class of the connection
con.id = length(n2.con)+1;
con.class = 'connection';
con.nClass = 3;

% obtain attributes from the source
con.source = n1.id;
con.sourceClass = n1.class;
con.nSourceClass = n1.nClass;   % numerical class (1: stimulus, 2: network)
con.sourceAxis = n1.f;
con.sourceAxisScale = n1.fspac;
con.nSourceAxisScale = n1.nFspac;
con.sourceAxisTick = n1.tick;
con.sourceN = n1.N;

% obtain attributes from the target
con.target = n2.id;
con.targetClass = n2.class;
con.nTargetClass = n2.nClass;
con.targetAxis = n2.f;
con.targetAxisScale = n2.fspac;
con.nTargetAxisScale = n2.nFspac;
con.targetAxisTick = n2.tick;
con.targetN = n2.N;

% parameters for the graphical display of the connection over time
con.type = '1freq';
con.dStep = 0;
con.sStep = 0;
con.no11 = 0;
con.tol  = .01;         % tolerance for fareyratio (used for 2freq)
con.phaseDisp = 0;
con.nType = 1; % numerical connection type
con.scale = 0;

% Connection parameters
F2 = repmat(n2.f, 1, n1.N);
F1 = repmat(n1.f', n2.N, 1);
F = (F1 + F2)/2;
con.F = ones(size(F));
con.w = 1;
con.lambda = 0;
con.mu1 = 0;
con.mu2 = 0;
con.kappa = 0;
con.e = 0;
con.learn = 0;
c21 = 197960;
bm2oc = diag(c21 * n2.f);
con.C0 = bm2oc;
con.C = con.C0;
n2.con{con.id} = con;

%% Set up the full network ================================================

% functions to calculate derivatives of networks and connections
model.zfun = @zdot;
model.cfun = @cdot;

% initialize model parameters
stimListAll = []; % list of all stimulus id's
stimList = []; % list of stimulus id's, only those used as source
netList = []; % list of network id's

% add stimulus to the model
sid = s.id;
s.z = s.x(:,1);   % initialize current state z
model.s{sid} = s;
stimListAll = [stimListAll sid];

% DSP parameters of the model
model.fs           = s.fs;
model.dt           = s.dt;
model.t            = s.t;
model.iSpan        = [1 length(s.t)];
model.iStep        = 1;

% add network 1 to the model
nid = n1.id;
model.n{nid} = n1;
netList = [netList nid];
Nt = ceil(length(s.t)/n1.sStep);
model.n{nid}.t = s.t(1:n1.sStep:length(s.t));
model.n{nid}.Z = zeros(length(n1.z), Nt);
model.n{nid}.Z(:,1) = n1.z0;

% add network 2 to the model
nid = n2.id;
model.n{nid} = n2;
netList = [netList nid];
Nt = ceil(length(s.t)/n2.sStep);
model.n{nid}.t = s.t(1:n2.sStep:length(s.t));
model.n{nid}.Z = zeros(length(n2.z), Nt);
model.n{nid}.Z(:,1) = n2.z0;

for j = netList
    con = model.n{j}.con{1};
    stimList = [stimList con.source];
end

model.stimList = stimList;
model.netList = netList;

M = model;

%% Run the full model =====================================================

zfun = M.zfun;
cfun = M.cfun;
iSpan = M.iSpan;
iStep = M.iStep;
h = M.dt;                   % step size
stimList = M.stimList;
netList = M.netList;

% Display stimulus and initial conditions if dStep > 0
stimulusLiveDisplay(M, 1, 0, M.t(1));

for nx = netList
    if M.n{nx}.dStep        
        networkLiveDisplay(M, nx, 0, M.t(1));
    end 
end

% Integration loop
for ix = iSpan(1) : iStep : iSpan(2)-iStep
    ind = ix; % time step for which to calculate k1     
    
    % Get Runge-Kutta k-values
    for kx = 1:4
        
        %% First update stimulus values
        if kx == 2 || kx == 4
            for sx = stimList
                M.s{sx}.z = stimulusRun(M.s{sx}, ind);
            end
        end
        
        %% Get k-values for each network
        for nx = netList
            M.n{nx}.k{kx} = h*zfun(M, nx);
            
            %% ... and for each learned connection to the network
            for cx = M.n{nx}.learnList
                M.n{nx}.con{cx}.k{kx} = h*cfun(M, nx, cx);
            end
        end
        
        %% Update z, C and ind for the next k-step
        switch kx
            case 1
                for nx = netList
                    M.n{nx}.zPrev = M.n{nx}.z;
                    M.n{nx}.z = M.n{nx}.zPrev + M.n{nx}.k{1}/2;
                    for cx = M.n{nx}.learnList
                        M.n{nx}.con{cx}.CPrev = M.n{nx}.con{cx}.C;
                        M.n{nx}.con{cx}.C = M.n{nx}.con{cx}.CPrev + M.n{nx}.con{cx}.k{1}/2;
                    end
                end
                ind = ix + iStep/2; % time step for k2 and k3
            case 2
                for nx = netList
                    M.n{nx}.z = M.n{nx}.zPrev + M.n{nx}.k{2}/2;
                    for cx = M.n{nx}.learnList
                        M.n{nx}.con{cx}.C = M.n{nx}.con{cx}.CPrev + M.n{nx}.con{cx}.k{2}/2;
                    end
                end
            case 3
                for nx = netList
                    M.n{nx}.z = M.n{nx}.zPrev + M.n{nx}.k{3};
                    for cx = M.n{nx}.learnList
                        M.n{nx}.con{cx}.C = M.n{nx}.con{cx}.CPrev + M.n{nx}.con{cx}.k{3};
                    end
                end
                ind = ix + iStep; % time step for k4
            case 4
                for nx = netList
                    M.n{nx}.z = M.n{nx}.zPrev + ...
                        (M.n{nx}.k{1} + 2*M.n{nx}.k{2} + 2*M.n{nx}.k{3} + M.n{nx}.k{4})/6;
                    if M.n{nx}.sStep && ~mod(ix, M.n{nx}.sStep)
                        M.n{nx}.Z(:,ix/M.n{nx}.sStep+1) = M.n{nx}.z;
                    end
                    if M.n{nx}.dStep && ~mod(ix, M.n{nx}.dStep)
                        networkLiveDisplay(M, nx, ix, M.t(ix));
                    end
                    for cx = M.n{nx}.learnList
                        M.n{nx}.con{cx}.C = M.n{nx}.con{cx}.CPrev + ...
                            (M.n{nx}.con{cx}.k{1} + 2*M.n{nx}.con{cx}.k{2} + 2*M.n{nx}.con{cx}.k{3} + M.n{nx}.con{cx}.k{4})/6;
                        if M.n{nx}.con{cx}.sStep && ~mod(ix, M.n{nx}.con{cx}.sStep)
                            M.n{nx}.con{cx}.C3(:,:,ix/M.n{nx}.con{cx}.sStep+1) = M.n{nx}.con{cx}.C;
                        end
                        if M.n{nx}.con{cx}.dStep && ~mod(ix, M.n{nx}.con{cx}.dStep)
                            connectionLiveDisplay(M, nx, cx, ix, M.t(ix));
                        end
                    end
                end
                
                for sx = stimList
                    if M.s{sx}.dStep && ~mod(ix, M.s{sx}.dStep)
                        stimulusLiveDisplay(M, sx, ix, M.t(ix));
                    end
                end
        end
    end
end

%%
figure(101);
plot(s.t, real(M.n{2}.Z(index,:)));zoom xon