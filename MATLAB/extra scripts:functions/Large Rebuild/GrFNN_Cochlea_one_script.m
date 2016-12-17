%% Decomposition of the Cochlear Model by Ed Large at al. (unpublished)
% this is the unidirectional version of the model

clear all
clc
close all

%% Define cochlear parameters =================================================
abm  = -412; aoc  =      0;
b1bm = 0;    b1oc =  -40816;
b2bm = 0;    b2oc =      0;
d1bm = 0;    d1oc =      0;
d2bm = 0;    d2oc =      0;
ebm  = 0;    eoc  =      0;

c21 = 197960;

%% Make a cochlea network =====================================================
% n1 = networkMake(1, 'hopf', abm, b1bm, b2bm, d1bm, d2bm, ebm, ...
%                     'log', 30, 10000, 831, ...
%                     'display', 50, 'save', 1, 'znaught', 0, 'noScale');

n1.id = 1;
n1.class = 'network';
n1.nClass = 2;

n1.con = {}; % cell array containing connections targeting this network
n1.learnList = []; % indices for learned connections (used in integrator)
models = {'hopf'};              % Can add to this array later
n1.model = [];                   % Initialize these to use isempty to error check later
n1.fspac = [];

n1.nFspac= 0;
n1.dStep = 0;                    % Initialize these to zero/empty in case not specified in varargin
n1.sStep = 0;
n1.tick  = [];
overrideInitialConditions = 0;
userScale = [];

alpha   = -412;
beta1   = 0;
beta2   = 0;
delta1  = 0;
delta2  = 0;
epsilon = 0;


n1.fspac = 'log';
lf   = 30;            % min freq
hf   = 10000;            % max freq
N    = 831;            % number of frequency steps
n1.N  = N;

n1.nFspac = 2;
n1.scale = 1;
n1.f  = logspace(log10(lf),log10(hf),N)';
n1.df = abs(log2(n1.f(2)/n1.f(1)));          % to scale for frequency density

userScale = 0;

n1.dStep = 50;

n1.sStep = 1;

overrideInitialConditions = 1;
z0 = 0;

n1.a  = alpha + 1i*2*pi.*n1.f;
n1.b1 = beta1 + 1i*delta1;
n1.b2 = beta2 + 1i*delta2;

n1.e = epsilon;

r0 = zeros(size(n1.a));
% r = spontAmp(real(n1.a(1)), real(n1.b1(1)), real(n1.b2(1)), n1.e);
% Find r* numerically
r = roots([n1.e*(real(n1.b2)-real(n1.b1)), 0, real(n1.b1)-real(n1.e)*real(n1.a(1)), 0, real(n1.a(1)), 0]);
r = real(unique(r(find(abs(imag(r)) < eps('single')))));
                    % only unique real values
r = r(find(r >=0)); % no negative amplitude

% Take only stable r*
drdotdr = real(n1.a(1)) + 3*n1.b1*r.^2 + (5*n1.e*n1.b2*r.^4-3*n1.e^2*n1.b2*r.^6)./((1-n1.e*r.^2).^2);
ind1 = find(drdotdr < 0);
ind2a = find(drdotdr == 0);
drdotdr = real(n1.a(1)) + 3*n1.b1*(r-eps('single')).^2 + (5*n1.e*n1.b2*(r-eps('single')).^4-3*n1.e^2*n1.b2*(r-eps('single')).^6)./((1-n1.e*(r-eps('single')).^2).^2);
ind2b = find(drdotdr < 0);
drdotdr = real(n1.a(1)) + 3*n1.b1*(r+eps('single')).^2 + (5*n1.e*n1.b2*(r+eps('single')).^4-3*n1.e^2*n1.b2*(r+eps('single')).^6)./((1-n1.e*(r+eps('single')).^2).^2);
ind2c = find(drdotdr < 0);
ind2 = intersect(ind2a,intersect(ind2b,ind2c));
r = r([ind1; ind2]);
r = sort(r,'descend');

% % ========================================================
% function drdotdr = slope(r, a, b1, b2, e)
% drdotdr = a + 3*b1*r.^2 + (5*e*b2*r.^4-3*e^2*b2*r.^6)./((1-e*r.^2).^2);

r0 = r(end)+r0;
r0 = r0+.01*randn(size(r0));
phi0 = 2*pi*randn(size(r0));
n1.z0 = r0.*exp(1i*2*pi*phi0);

n1.z0 = repmat(z0,N,1);

n1.z = n1.z0;

% n2 = networkMake(2, 'hopf', aoc, b1oc, b2oc, d1oc, d2oc, eoc, ...
%                     'log', 30, 10000, 831, ...
%                     'display', 50, 'save', 1, 'znaught', 0, 'noScale');

n2.id = 2;
n2.class = 'network';
n2.nClass = 2;

n2.con = {}; % cell array containing connections targeting this network
n2.learnList = []; % indices for learned connections (used in integrator)
models = {'hopf'};              % Can add to this array later
n2.model = [];                   % Initialize these to use isempty to error check later
n2.fspac = [];

n2.nFspac= 0;
n2.dStep = 0;                    % Initialize these to zero/empty in case not specified in varargin
n2.sStep = 0;
n2.tick  = [];
overrideInitialConditions = 0;
userScale = [];

alpha   = 0;
beta1   = -40816;
beta2   = 0;
delta1  = 0;
delta2  = 0;
epsilon = 0;

n2.fspac = 'log';
lf   = 30;            % min freq
hf   = 10000;            % max freq
N    = 831;            % number of frequency steps
n2.N  = N;

n2.nFspac = 2;
n2.scale = 1;
n2.f  = logspace(log10(lf),log10(hf),N)';
n2.df = abs(log2(n2.f(2)/n2.f(1)));          % to scale for frequency density

userScale = 0;

n2.dStep = 50;

n2.sStep = 1;

overrideInitialConditions = 1;
z0 = 0;

n2.a  = alpha + 1i*2*pi.*n2.f;
n2.b1 = beta1 + 1i*delta1;
n2.b2 = beta2 + 1i*delta2;

n2.e = epsilon;

r0 = zeros(size(n2.a));
r = spontAmp(real(n2.a(1)), real(n2.b1(1)), real(n2.b2(1)), n2.e);
r0 = r(end)+r0;
r0 = r0+.01*randn(size(r0));
phi0 = 2*pi*randn(size(r0));
n2.z0 = r0.*exp(1i*2*pi*phi0);

n2.z0 = repmat(z0,N,1);

n2.z = n2.z0;

%% Make a stimulus ============================================================
F = dB2Pa(10);
freqs = 1850;
% index = freqToIndex(n1, freqs);
[NF, FREQ] = meshgrid(n1.f, freqs);
[~, index] = min(abs(log2(NF ./ FREQ)), [], 2);

freqs = n1.f(index);

% s = stimulusMake(1, 'fcn', [0 .10], 100000, {'exp'}, freqs, F, 0, ...
%                  'ramp', 0.010, 1, 'display', 200);

id = 1;
type = 'fcn';

s.id = id;
s.class = 'stimulus';
s.nClass = 1; % numerical class
s.dStep = 200;
s.dispChan = 1;
s.f = [];
s.fspac = '';
s.nFspac = 0;
s.tick = [];

s.type = type; 

s.ts = [0 .10];         % time spans
s.fs = 100000;         % sampling frequency
s.dt = 1/s.fs;
s.carrier = {'exp'};    % carrier waveforms
s.fc = freqs;         % carrier frequencies
s.ac = F;         % carrier amplitudes
s.sc = .02;                 % ramp time in sec of beginning and end of stim timecourses
s.sp = 1;                   % ramp strength exponent, ie, larger than 1, sudden onset, smaller than one, gradual

varargin = cell([[0 .10], 100000, {'exp'}, freqs, F, 0, ...
                  'ramp', 0.010, 1, 'display', 200]);

s = stimulusParser(s, varargin{:});

s.t  = min(min(s.ts)):s.dt:max(max(s.ts));
s.x  = zeros(size(s.t));

for a = 1:size(s.ts,1)                % For each section of the signal
    t0n = find(s.t <= s.ts(a,1), 1, 'last' );
    tfn = find(s.t <= s.ts(a,2), 1, 'last' );
    %     t  = s.t(t0n:tfn);
    t  = s.t(1:tfn-t0n+1);
    
    temp = stimulusFcn(t, s, a);  %quicker to write temp vector than keep indexing into s.x
    
    if any(s.sc) && ~strcmp(s.carrier{1}, 'pls') %kludgy but works for now
        temp = stimulusRamp(temp, s.sc(a), s.sp(a), s.fs);
    end
    
    if isfield(s, 'filtstim') && ~isempty(s.filtstim{a,1}) && ~isempty(s.filtstim{a,2})
        temp = filter(s.filtstim{a,1},s.filtstim{a,2},temp);
    end
    
    if isfield(s, 'mask')
        noise = rand(size(temp))*2-1;
        if isfield(s, 'filtmask') && ~isempty(s.filtmask{a,1}) && ~isempty(s.filtmask{a,2})
            noise = filter(s.filtmask{a,1},s.filtmask{a,2},noise);
        end
        noise = noise * (1/rootMeanSquare(noise)) * 10^(-s.mask(a)/20) * rootMeanSquare(temp);    %interpret noise input as SNR in dB with reference to stim
        
        noise = stimulusRamp(noise, s.sc(a), s.sp(a), s.fs);
        temp = temp + noise;
    end
    
    s.x(t0n:tfn) = s.x(t0n:tfn) + temp;            %write temp vector to approp. part of s.x
    
end

if isfield(s, 'maskall')
    noise = rand(size(s.x))*2-1;
    if isfield(s, 'filtmaskall') && ~isempty(s.filtmaskall{1}) && ~isempty(s.filtmaskall{2})
        noise = filter(s.filtmaskall{1},s.filtmaskall{2},noise);
    end
    noise = noise * (1/rootMeanSquare(noise)) * 10^(-s.maskall/20) * rootMeanSquare(s.x);
    
    noise = stimulusRamp(noise, s.sc(1), s.sp(1), s.fs);
    s.x = s.x + noise;
end

s.lenx = size(s.x, 2);	% stimulus length
s.N = size(s.x, 1);     % number of stimulus channels

Fs = 100e3;
tdres = 1/Fs;
fp    = 1e3;  % prewarping frequency 1 kHz 
TWOPI = 2*pi;
C     = TWOPI*fp/tan((TWOPI/2)*fp*tdres);

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

% s.x = midearfilt(s.x, s.fs);

% n1 = connectAdd(s, n1, 1, 'noScale');
con.id = length(n1.con)+1;
con.class = 'connection';
con.nClass = 3;

con.source = s.id;
con.sourceClass = s.class;
con.nSourceClass = s.nClass;   % numerical class (1: stimulus, 2: network)
con.sourceAxis = s.f;
con.sourceAxisScale = s.fspac;
con.nSourceAxisScale = s.nFspac;
con.sourceAxisTick = s.tick;
con.sourceN = s.N;

con.target = n1.id;
con.targetClass = n1.class;
con.nTargetClass = n1.nClass;
con.targetAxis = n1.f;
con.targetAxisScale = n1.fspac;
con.nTargetAxisScale = n1.nFspac;
con.targetAxisTick = n1.tick;
con.targetN = n1.N;

con.type = '1freq';

con.dStep = 0;
con.sStep = 0;
con.no11 = 0;
con.tol  = .01;         % tolerance for fareyratio (used for 2freq)
con.phaseDisp = 0;
w = 1;           % Initialize these in case not specified in varargin
learn = 0;
lambda = 0;
mu1 = 0;
mu2 = 0;
epsilon = 0;
kappa = 0;
userScale = 0;

F2 = repmat(n1.f, 1, s.N);
F1 = F2;

F = (F1 + F2)/2;
con.nType = 1; % numerical connection type

con.scale = userScale;

con.F = ones(size(F));
con.w = w;
con.lambda = lambda;
con.mu1 = mu1;
con.mu2 = mu2;
con.kappa = kappa;

con.e = epsilon;
con.learn = learn;
con.C0 = ones(size(F));
con.C = con.C0;

mask = 1;

con.lambda = con.lambda .* mask;
con.mu1 = con.mu1 .* mask;
con.mu2 = con.mu2 .* mask;
con.kappa = con.kappa .* mask;
con.C0    = con.C0    .* mask;
con.C     = con.C     .* mask;
n1 = n1;
n1.con{con.id} = con;

%% Add connections from bm to oc
bm2oc = diag(c21 * n2.f);

% n2    = connectAdd(n1, n2, bm2oc, 'type', '1freq', 'noScale');

con.id = length(n2.con)+1;
con.class = 'connection';
con.nClass = 3;

con.source = n1.id;
con.sourceClass = n1.class;
con.nSourceClass = n1.nClass;   % numerical class (1: stimulus, 2: network)
con.sourceAxis = n1.f;
con.sourceAxisScale = n1.fspac;
con.nSourceAxisScale = n1.nFspac;
con.sourceAxisTick = n1.tick;
con.sourceN = n1.N;

con.target = n2.id;
con.targetClass = n2.class;
con.nTargetClass = n2.nClass;
con.targetAxis = n2.f;
con.targetAxisScale = n2.fspac;
con.nTargetAxisScale = n2.nFspac;
con.targetAxisTick = n2.tick;
con.targetN = n2.N;

con.type = 'Allfreq'; 

con.dStep = 0;
con.sStep = 0;
con.no11 = 0;
con.tol  = .01;         % tolerance for fareyratio (used for 2freq)
con.phaseDisp = 0;
w = 1;           % Initialize these in case not specified in varargin
learn = 0;
lambda = 0;
mu1 = 0;
mu2 = 0;
epsilon = 0;
kappa = 0;
userScale = [];

con.type = '1freq';

userScale = 0;

F2 = repmat(n2.f, 1, n1.N);
F1 = repmat(n1.f', n2.N, 1);

F = (F1 + F2)/2;
con.nType = 1; % numerical connection type

con.scale = userScale; 

con.F = ones(size(F));
con.w = w;
con.lambda = lambda;
con.mu1 = mu1;
con.mu2 = mu2;
con.kappa = kappa;

con.e = epsilon;
con.learn = learn;

con.C0 = bm2oc;
con.C = con.C0;

mask = 1;

con.lambda = con.lambda .* mask;
con.mu1 = con.mu1 .* mask;
con.mu2 = con.mu2 .* mask;
con.kappa = con.kappa .* mask;
con.C0    = con.C0    .* mask;
con.C     = con.C     .* mask;

n2 = n2;
n2.con{con.id} = con;

%% Run the network
% M = modelMake(@zdot, @cdot, s, n1, n2);

model.zfun = @zdot;
model.cfun = @cdot;
ind = 1;

model.odefun = @odeRK4fs;   % default ode function
model.usingGpu = 0;

stimListAll = []; % list of all stimulus id's
stimList = []; % list of stimulus id's, only those used as source
netList = []; % list of network id's
fs = [];
dt = [];
ts = [];
t = [];

sid = s.id;
s.N = size(s.x, 1);   % in case extra channel was added
s.z = s.x(:,1);   % initialize current state z
model.s{sid} = s;
stimListAll = [stimListAll sid];

fs = s.fs;
dt = s.dt;
ts = s.ts;
t = s.t;

model.fs           = fs;
model.dt           = dt;
model.t            = t;
model.iSpan        = [1 length(t)];
model.iStep        = 1;

nid = n1.id;
model.n{nid} = n1;
netList = [netList nid];

Nt = ceil(length(t)/n1.sStep);
model.n{nid}.t = t(1:n1.sStep:length(t));
model.n{nid}.Z = zeros(length(n1.z), Nt);
model.n{nid}.Z(:,1) = n1.z0;

nid = n2.id;
model.n{nid} = n2;
netList = [netList nid];

Nt = ceil(length(t)/n2.sStep);
model.n{nid}.t = t(1:n2.sStep:length(t));
model.n{nid}.Z = zeros(length(n2.z), Nt);
model.n{nid}.Z(:,1) = n2.z0;

netList = sort(netList);
stimListAll = sort(stimListAll);

for j = netList
    for k = 1:length(model.n{j}.con)
        con = model.n{j}.con{k};
        if con.nSourceClass == 1
            if ~ismember(con.source, stimListAll)
                error(['Input to Network ' num2str(j) ' is missing (Stimulus ' num2str(k) ')'])
            else
                stimList = [stimList con.source];
            end
        end
        if con.nSourceClass == 2
            if ~ismember(con.source, netList)
                error(['Input to Network ' num2str(j) ' is missing (Network ' num2str(k) ')'])
            end
        end
    end
end

model.stimListAll = stimListAll;
model.stimList = sort(stimList);
model.netList = netList;

model.iSpan    = single(model.iSpan);
model.iStep    = single(model.iStep);
model.dt       = single(model.dt);
for sx = model.stimList
    model.s{sx}.x = castCS(model.s{sx}.x);
end

for nx = model.netList
    model.n{nx}.z0 = castCS(model.n{nx}.z0);
    model.n{nx}.z  = castCS(model.n{nx}.z);
    model.n{nx}.Z  = castCS(model.n{nx}.Z);
    model.n{nx}.a  = castCS(model.n{nx}.a);
    model.n{nx}.b1 = castCS(model.n{nx}.b1);
    model.n{nx}.b2 = castCS(model.n{nx}.b2);
    model.n{nx}.e  = castCS(model.n{nx}.e);
    for cx = 1:length(model.n{nx}.con)
        if model.n{nx}.con{cx}.learn
            model.n{nx}.con{cx}.C3 = castCS(model.n{nx}.con{cx}.C3);
		end
		% DMR Do not castCS on NUM, NUM1, NUM2, and DEN. This will crash the C code (mex).
		% But we _do_ want singles.
        if isfield(model.n{nx}.con{cx},'NUM')
            model.n{nx}.con{cx}.NUM = single(model.n{nx}.con{cx}.NUM);
            model.n{nx}.con{cx}.DEN = single(model.n{nx}.con{cx}.DEN);
        end
        if isfield(model.n{nx}.con{cx},'NUM1')
            model.n{nx}.con{cx}.NUM1 = single(model.n{nx}.con{cx}.NUM1);
            model.n{nx}.con{cx}.NUM2 = single(model.n{nx}.con{cx}.NUM2);
            model.n{nx}.con{cx}.DEN  = single(model.n{nx}.con{cx}.DEN);
        end
        model.n{nx}.con{cx}.C      = castCS(model.n{nx}.con{cx}.C);
        model.n{nx}.con{cx}.C0     = castCS(model.n{nx}.con{cx}.C0);
        model.n{nx}.con{cx}.w      = castCS(model.n{nx}.con{cx}.w);
        model.n{nx}.con{cx}.lambda = castCS(model.n{nx}.con{cx}.lambda);
        model.n{nx}.con{cx}.mu1    = castCS(model.n{nx}.con{cx}.mu1);
        model.n{nx}.con{cx}.mu2    = castCS(model.n{nx}.con{cx}.mu2);
        model.n{nx}.con{cx}.kappa  = castCS(model.n{nx}.con{cx}.kappa);
        model.n{nx}.con{cx}.e      = castCS(model.n{nx}.con{cx}.e);
    end
end

M = model;

tic;
M = odeRK4fs(M);
clc; toc;

%%
figure(101);
plot(s.t, real(M.n{2}.Z(index,:)));zoom xon                