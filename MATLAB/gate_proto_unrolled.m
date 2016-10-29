close all
clear all
clc

% Include the GrFNN toolbox
addpath(genpath('~/Documents/GrFNNToolbox/'))

%% Network parameters
% alpha = 1; beta1 = -1; beta2 = -1000; neps = 1; % Limit Cycle
alpha = 1; beta1 = -1; beta2 = -0.1; delta1 = 0; delta2 = 0; neps = 1; % Limit Cycle

%% Parameter sets for Hebbian plasiticity
w = .05;
lambda =  -.1; mu1 =  0; mu2 =  0; ceps =  4; kappa = 1; % Linear learning rule
% lambda =   0; mu1 = -1; mu2 = -50; ceps =  4; kappa = 1; % Critical learning rule
% lambda =   0; mu1 = -1; mu2 = -50; ceps = 16; kappa = 1; % Critical, stronger nonlinearity
% lambda = .001; mu1 = -1; mu2 = -50; ceps = 16; kappa = 1; % Supercritical learning rule

clear all
close all
clc

% Include the GrFNN toolbox
addpath(genpath('~/Documents/GrFNNToolbox/'))

%% Network parameters
% alpha = 1; beta1 = -1; beta2 = -1000; neps = 1; % Limit Cycle
alpha = 1; beta1 = -1; beta2 = -0.1; delta1 = 0; delta2 = 0; neps = 1; % Limit Cycle

%% Parameter sets for Hebbian plasiticity
w = .05;
lambda =  -.1; mu1 =  0; mu2 =  0; ceps =  4; kappa = 1; % Linear learning rule
% lambda =   0; mu1 = -1; mu2 = -50; ceps =  4; kappa = 1; % Critical learning rule
% lambda =   0; mu1 = -1; mu2 = -50; ceps = 16; kappa = 1; % Critical, stronger nonlinearity
% lambda = .001; mu1 = -1; mu2 = -50; ceps = 16; kappa = 1; % Supercritical learning rule

%% Make the model
s = stimulusMake(1, 'fcn', [0 1], 4000, {'exp'}, 1, 0);

s.x = zeros(size(s.x));
for i=1:500:length(s.x)
    s.x(i:(i+100)) = 10;
end
% plot(s.x)
% pause

%%
n.id = 1;
n.class = 'network';
n.nClass = 2; % numerical class

n.con = {}; % cell array containing connections targeting this network
n.learnList = []; % indices for learned connections (used in integrator)

models = {'hopf'};              % Can add to this array later
n.model = [];                   % Initialize these to use isempty to error check later
n.fspac = [];

n.nFspac= 0;
n.dStep = 0;                    % Initialize these to zero/empty in case not specified in varargin
n.sStep = 0;
overrideInitialConditions = 0;
n.tick  = [];

varargin = {1, 'hopf', alpha, beta1,  beta2, 0, 0, neps, ...
    'lin', 15, 19, 20, 'save', 1, ...
    'display', 10};

for i = 1:length(varargin)
    
    if any(strcmpi(varargin{i},models)) && length(varargin) > i + 5 && ~ischar(varargin{i+1}) && ~ischar(varargin{i+2}) && ~ischar(varargin{i+3}) && ~ischar(varargin{i+4}) && ~ischar(varargin{i+5}) && ~ischar(varargin{i+6})
        
        n.model = lower(varargin{i});
        
        switch varargin{i}
            
            case {'hopf'} % 'hopp'? Let's make it 'hopf' ...
                alpha   = varargin{i+1};
                beta1   = varargin{i+2};
                beta2   = varargin{i+3};
                delta1  = varargin{i+4};
                delta2  = varargin{i+5};
                epsilon = varargin{i+6};
                
        end
    end
    
    if ischar(varargin{i}) && any(strcmpi(varargin{i}(1:3),{'lin' 'log'})) && length(varargin) > i + 2 && isscalar(varargin{i+1}) && isscalar(varargin{i+2}) && isscalar(varargin{i+3})
        
        n.fspac = lower(varargin{i}(1:3));
        % Assign nFspac integer value based on fspac string value
        if strcmpi(n.fspac, 'lin')
            n.nFspac = 1;
        else
            n.nFspac = 2;
        end
        
        lf   = varargin{i+1};            % min freq
        hf   = varargin{i+2};            % max freq
        N    = varargin{i+3};            % number of frequency steps
        n.N  = N;
        switch n.nFspac
            
            case 1 % linear spacing
                n.f  = linspace(lf,hf,N)';
                if N > 1
                    n.df = abs(n.f(2)-n.f(1)); % to scale for frequency density
                else
                    n.df = 1;
                end
                
            case 2 % log spacing
                n.f  = logspace(log10(lf),log10(hf),N)';
                if N > 1
                    n.df = abs(log2(n.f(2)/n.f(1)));          % to scale for frequency density
                else
                    n.df = 1;
                end
        end
        n.f   = single(n.f);
    end
     
    if ischar(varargin{i}) && strcmpi(varargin{i}(1:3),'dis') && length(varargin) > i && isscalar(varargin{i+1})
        
        n.dStep = varargin{i+1};
        
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i}(1:3),'sav') && length(varargin) > i && isscalar(varargin{i+1})
        
        n.sStep = varargin{i+1};
        
    end
    
end

switch n.nFspac
    
    case 1 % linear spacing
        n.a  = single(alpha + 1i*2*pi.*n.f);
        n.b1 = single(beta1 + 1i*delta1);
        n.b2 = single(beta2 + 1i*delta2);
        
    case 2 % log spacing
        n.a  = single(alpha + 1i*2*pi  ).*n.f;  % Redefinition of a, b1 & b2
        n.b1 = single(beta1 + 1i*delta1).*n.f;
        n.b2 = single(beta2 + 1i*delta2).*n.f;
end

n.e = single(epsilon);

switch n.model
    case {'hopf', 'hopp'}
        r0 = zeros(size(n.a));
        r = spontAmp(real(n.a(1)), real(n.b1(1)), real(n.b2(1)), n.e);
        r0 = r(end)+r0;
        r0 = r0+.01*randn(size(r0));
        phi0 = 2*pi*randn(size(r0));
        n.z0 = r0.*exp(1i*2*pi*phi0);
        
    otherwise
        error('Unknown model type: %s', n.model);
        
end

if overrideInitialConditions
    if numel(z0) == 1
        n.z0 = repmat(z0,N,1);
    elseif all(size(z0) == [1 N])
        n.z0 = z0.';
    elseif all(size(z0) == [N 1])
        n.z0 = z0;
    else
        error('Length of initial conditions must be 1 or N')
    end        
end

n.z = n.z0;

n1 = n;
n2 = n;
C = [];

varargin = {'weight', w, 'type', 'all2freq', ...
    'learn', lambda, mu1, mu2, ceps, kappa, ...
    'display', 10,'phasedisp', 'save', 500};

% n1 is the source (stimulus or network), n2 is the target network
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

%% Initialize parameters

if con.nSourceClass == 1
    con.type = '1freq';     % default connection type for stimulus source
else
    con.type = 'Allfreq';   % default connection type for network source
end
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

%% Parse input
types = {'1freq' '2freq' '3freq' '3freqAll' 'All2freq' 'Allfreq' 'active'};         % Can add to this array later; default to 'Allfreq' for now

for i = 1:length(varargin)
    
    if strcmpi(varargin{i},'learn') && length(varargin) > i + 4 && ~ischar(varargin{i+1}) && ~ischar(varargin{i+2}) && ~ischar(varargin{i+3}) && ~ischar(varargin{i+4}) && ~ischar(varargin{i+5})
        learn = 1;
        lambda   = varargin{i+1};
        mu1      = varargin{i+2};
        mu2      = varargin{i+3};
        epsilon  = varargin{i+4};
        kappa    = varargin{i+5};
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i}(1:3),'typ') && length(varargin) > i && ischar(varargin{i+1}) && any(strcmpi(varargin{i+1},types))
        con.type = varargin{i+1};
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i}(1:3),'2fr') && length(varargin) > i && isscalar(varargin{i+1})
        con.tol = varargin{i+1};
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i}(1:3),'wei') && length(varargin) > i && isnumeric(varargin{i+1})
        w = varargin{i+1};
    end
    
    
    if ischar(varargin{i}) && strcmpi(varargin{i}(1:3),'dis') && length(varargin) > i && isscalar(varargin{i+1})
        con.dStep = varargin{i+1};
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i}(1:3),'sav') && length(varargin) > i && isscalar(varargin{i+1})
        con.sStep = varargin{i+1};
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i},'no11')
        con.no11 = 1;
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i},'phasedisp')
        con.phaseDisp = 1;
    end
    
    if ischar(varargin{i}) && ~strcmpi(varargin{i},'learn') && ~strcmpi(varargin{i}(1:3),'typ') && ~strcmpi(varargin{i}(1:3),'wei') && ~strcmpi(varargin{i}(1:3),'dis') && ~strcmpi(varargin{i}(1:3),'sav') && ~any(strcmpi(varargin{i},types)) && ~strcmpi(varargin{i},'no11') && ~strcmpi(varargin{i},'phasedisp')
        error(['Unrecognized input to connectAdd: ' varargin{i}])
    end
    
end

%% Connection types

F2 = repmat(n2.f, 1, n1.N);

if n1.nClass == 1 % if stimulus source
    F1 = F2;
    if isempty(n1.f) && ismember(lower(con.type), {'2freq', '3freq', '3freqall'})
        error(['Connection type ''' lower(con.type) ''' not available for stimulus source with no frequency gradient'])
    end
else
    F1 = repmat(n1.f', n2.N, 1);
end

switch lower(con.type)
    
    case '1freq' % one-frequency monomials
        F = (F1 + F2)/2;
        con.nType = 1; % numerical connection type
        
    case '2freq' % two-frequency monomials
        R = F2./F1;
        sz = size(R);
        [NUM, DEN] = fareyratio(R(:)', con.tol);
        NUM = reshape(NUM,sz);
        DEN = reshape(DEN,sz);
        con.NUM = NUM;
        con.DEN = DEN;
        F = (NUM.*F1 + DEN.*F2)./(NUM + DEN);
        if n1.nClass == 2 % if network source
            con.epsn = n1.e.^((NUM+DEN-2)/2); % used in zdot
        end
        con.epsc = epsilon.^((NUM+DEN-2)/2);  % used in cdot
        con.nType = 2;
        
    case '3freq' % three-frequency monomials ALL MONONMIALS, UP TO SPECIFIED ORDER
        if n1.id == n2.id && n1.nClass == 2     % if internal connection
            [X1i, X2i, Zi, NUM1, NUM2, DEN, CON1, CON2, mask] = threeFreqMatsAll(n1.f);
        else
            [X1i, X2i, Zi, NUM1, NUM2, DEN, CON1, CON2, mask] = threeFreqMatsAll(n1.f, n2.f);
        end
        con.IDX1 = X1i;
        con.IDX2 = X2i;
        con.IDXZ  = Zi;
        con.CON1 = CON1;
        con.CON2 = CON2;
        con.NUM1 = NUM1;
        con.NUM2 = NUM2;
        con.DEN = DEN;
        con.mask = mask;
        if n1.nClass == 1 % if stimulus source
            F = repmat(n2.f, 1, size(DEN, 2));
        else
            F = (abs(NUM1).*n1.f(X1i) + abs(NUM2).*n1.f(X2i) + DEN.*n2.f(Zi)) ...
                ./(abs(NUM1) + abs(NUM2) + DEN);
            con.epsn = n1.e.^((NUM1+NUM2+DEN-2)/2); % used in zdot
        end
        con.epsc = epsilon.^((NUM1+NUM2+DEN-2)/2);  % used in cdot
        con.nType = 3;
        
    case '3freqall' % three-frequency monomials ALL FREQUENCIES, LOWEST ORDER MONOMIAL
        if n1.id == n2.id && n1.nClass == 2     % if internal connection
            [X1i, X2i, Zi] = inputShapeInternal(n1.N);
            [NUM1,  NUM2,  DEN] = inputExponentsInternal(n1.f, X1i, X2i, Zi);
        else
            [X1i, X2i, Zi] = inputShapeOther(n1.N, n2.N);
            [NUM1,  NUM2,  DEN] = inputExponentsOther(n1.f, n2.f, X1i, X2i, Zi);
        end
        con.IDX1 = X1i;
        con.IDX2 = X2i;
        con.IDXZ  = Zi;
        con.CON1 = (NUM1<0);
        con.CON2 = (NUM2<0);
        con.NUM1 = abs(NUM1);
        con.NUM2 = abs(NUM2);
        con.DEN = DEN;
        if n1.nClass == 1 % if stimulus source
            F = repmat(n2.f, 1, size(DEN, 2));
        else
            F = (abs(NUM1).*n1.f(X1i) + abs(NUM2).*n1.f(X2i) + DEN.*n2.f(Zi)) ...
                ./(abs(NUM1) + abs(NUM2) + DEN);
            con.epsn = n1.e.^((NUM1+NUM2+DEN-2)/2); % used in zdot
        end
        con.epsc = epsilon.^((NUM1+NUM2+DEN-2)/2);  % used in cdot
        con.nType = 4;
        
    case 'active' % Full series of active nonlinearities
        F = (2*F1.*F2)./(F1+F2);
        con.nType = 5;
        
    case 'all2freq' % Full series of resonant monomials
        F = (2*F1.*F2)./(F1+F2);
        con.nType = 6;
        
    case 'allfreq' % Full series of resonant monomials including conjugates
        F = (2*F1.*F2)./(F1+F2);
        con.nType = 7;
        
end

%% Scaling factors

if n2.nFspac==2 % log spacing
    con.F = F;
    con.w = w.*n2.f;
else
    F = ones(size(F));
    con.w = w;
end

%% Scale learning parameters

% if any(kappa) % JCK: replaces if sum(sum(kappa))
%     con.learn = 1;
% else
%     con.learn = 0;
% end

con.learn = learn;

if n2.nFspac==2 % log spacing
    con.lambda = lambda.*F;
    con.mu1 = mu1.*F;
    con.mu2 = mu2.*F;
    con.kappa = kappa.*F;
else
    con.lambda = lambda;
    con.mu1 = mu1;
    con.mu2 = mu2;
    con.kappa = kappa;
end
con.e = epsilon;

%% Connection State and Initial Conditions
%       C0: initial conditions
%        C: connection matrix (instantaneous)
%        t: times saved in 3D memory matrix
%       C3: state memory (3D matrix: frequency x frequency x time)

if isempty(C)
    if ~con.learn
        error('Connection values must be specified as C (3rd input argument) if not learning')
    end
    A0 = zeros(size(F));
    A = spontAmp(real(con.lambda(1,1)), real(con.mu1(1,1)), real(con.mu2(1,1)), con.e);
    A0 = A0+min(A);
    A0 = A0.*(1 +.01*randn(size(A0)));
    theta0 = randn(size(A0));
    con.C0 = A0.*exp(1i*2*pi*theta0);
else
    if isscalar(C)
        con.C0 = C*ones(size(F));
    elseif ~isequal(size(C), size(F))
        error('Incorrect size of C (3rd input argument)')
    else
        con.C0 = C;
    end
end

con.C = con.C0;

%% Masking
%      This process operates on the initial condition and parameter
%      matrices to hange the way things work (i.e.
%
%      mask = ~eye(size(C));          % Don't learn self-connections
%                                       This could be inverted guassian
%                                       too

mask = 1;
if strcmpi(con.type, '2freq') && con.no11
    mask = ~(con.NUM == 1 & con.DEN == 1);
elseif strcmpi(con.type, '3freq')
    mask = con.mask;
elseif con.source == con.target && con.nSourceClass == 2
    mask = ~eye(size(con.C));          % ... Won't learn self-connections
end

con.lambda = con.lambda .* mask;
con.mu1 = con.mu1 .* mask;
con.mu2 = con.mu2 .* mask;
con.kappa = con.kappa .* mask;
con.C0    = con.C0    .* mask;
con.C     = con.C     .* mask;

%% Return network n2
n = n2;
n.con{con.id} = con;

%% Add index to learnList if learn
if con.learn
    n.learnList = [n.learnList con.id];
end

n1 = s;
n2 = n;
C = 1;

varargin = {};

% n1 is the source (stimulus or network), n2 is the target network
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

%% Initialize parameters

if con.nSourceClass == 1
    con.type = '1freq';     % default connection type for stimulus source
else
    con.type = 'Allfreq';   % default connection type for network source
end
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

%% Parse input
types = {'1freq' '2freq' '3freq' '3freqAll' 'All2freq' 'Allfreq' 'active'};         % Can add to this array later; default to 'Allfreq' for now

for i = 1:length(varargin)
    
    if strcmpi(varargin{i},'learn') && length(varargin) > i + 4 && ~ischar(varargin{i+1}) && ~ischar(varargin{i+2}) && ~ischar(varargin{i+3}) && ~ischar(varargin{i+4}) && ~ischar(varargin{i+5})
        learn = 1;
        lambda   = varargin{i+1};
        mu1      = varargin{i+2};
        mu2      = varargin{i+3};
        epsilon  = varargin{i+4};
        kappa    = varargin{i+5};
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i}(1:3),'typ') && length(varargin) > i && ischar(varargin{i+1}) && any(strcmpi(varargin{i+1},types))
        con.type = varargin{i+1};
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i}(1:3),'2fr') && length(varargin) > i && isscalar(varargin{i+1})
        con.tol = varargin{i+1};
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i}(1:3),'wei') && length(varargin) > i && isnumeric(varargin{i+1})
        w = varargin{i+1};
    end
    
    
    if ischar(varargin{i}) && strcmpi(varargin{i}(1:3),'dis') && length(varargin) > i && isscalar(varargin{i+1})
        con.dStep = varargin{i+1};
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i}(1:3),'sav') && length(varargin) > i && isscalar(varargin{i+1})
        con.sStep = varargin{i+1};
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i},'no11')
        con.no11 = 1;
    end
    
    if ischar(varargin{i}) && strcmpi(varargin{i},'phasedisp')
        con.phaseDisp = 1;
    end
    
    if ischar(varargin{i}) && ~strcmpi(varargin{i},'learn') && ~strcmpi(varargin{i}(1:3),'typ') && ~strcmpi(varargin{i}(1:3),'wei') && ~strcmpi(varargin{i}(1:3),'dis') && ~strcmpi(varargin{i}(1:3),'sav') && ~any(strcmpi(varargin{i},types)) && ~strcmpi(varargin{i},'no11') && ~strcmpi(varargin{i},'phasedisp')
        error(['Unrecognized input to connectAdd: ' varargin{i}])
    end
    
end

%% Connection types

F2 = repmat(n2.f, 1, n1.N);

if n1.nClass == 1 % if stimulus source
    F1 = F2;
    if isempty(n1.f) && ismember(lower(con.type), {'2freq', '3freq', '3freqall'})
        error(['Connection type ''' lower(con.type) ''' not available for stimulus source with no frequency gradient'])
    end
else
    F1 = repmat(n1.f', n2.N, 1);
end

switch lower(con.type)
    
    case '1freq' % one-frequency monomials
        F = (F1 + F2)/2;
        con.nType = 1; % numerical connection type
        
    case '2freq' % two-frequency monomials
        R = F2./F1;
        sz = size(R);
        [NUM, DEN] = fareyratio(R(:)', con.tol);
        NUM = reshape(NUM,sz);
        DEN = reshape(DEN,sz);
        con.NUM = NUM;
        con.DEN = DEN;
        F = (NUM.*F1 + DEN.*F2)./(NUM + DEN);
        if n1.nClass == 2 % if network source
            con.epsn = n1.e.^((NUM+DEN-2)/2); % used in zdot
        end
        con.epsc = epsilon.^((NUM+DEN-2)/2);  % used in cdot
        con.nType = 2;
        
    case '3freq' % three-frequency monomials ALL MONONMIALS, UP TO SPECIFIED ORDER
        if n1.id == n2.id && n1.nClass == 2     % if internal connection
            [X1i, X2i, Zi, NUM1, NUM2, DEN, CON1, CON2, mask] = threeFreqMatsAll(n1.f);
        else
            [X1i, X2i, Zi, NUM1, NUM2, DEN, CON1, CON2, mask] = threeFreqMatsAll(n1.f, n2.f);
        end
        con.IDX1 = X1i;
        con.IDX2 = X2i;
        con.IDXZ  = Zi;
        con.CON1 = CON1;
        con.CON2 = CON2;
        con.NUM1 = NUM1;
        con.NUM2 = NUM2;
        con.DEN = DEN;
        con.mask = mask;
        if n1.nClass == 1 % if stimulus source
            F = repmat(n2.f, 1, size(DEN, 2));
        else
            F = (abs(NUM1).*n1.f(X1i) + abs(NUM2).*n1.f(X2i) + DEN.*n2.f(Zi)) ...
                ./(abs(NUM1) + abs(NUM2) + DEN);
            con.epsn = n1.e.^((NUM1+NUM2+DEN-2)/2); % used in zdot
        end
        con.epsc = epsilon.^((NUM1+NUM2+DEN-2)/2);  % used in cdot
        con.nType = 3;
        
    case '3freqall' % three-frequency monomials ALL FREQUENCIES, LOWEST ORDER MONOMIAL
        if n1.id == n2.id && n1.nClass == 2     % if internal connection
            [X1i, X2i, Zi] = inputShapeInternal(n1.N);
            [NUM1,  NUM2,  DEN] = inputExponentsInternal(n1.f, X1i, X2i, Zi);
        else
            [X1i, X2i, Zi] = inputShapeOther(n1.N, n2.N);
            [NUM1,  NUM2,  DEN] = inputExponentsOther(n1.f, n2.f, X1i, X2i, Zi);
        end
        con.IDX1 = X1i;
        con.IDX2 = X2i;
        con.IDXZ  = Zi;
        con.CON1 = (NUM1<0);
        con.CON2 = (NUM2<0);
        con.NUM1 = abs(NUM1);
        con.NUM2 = abs(NUM2);
        con.DEN = DEN;
        if n1.nClass == 1 % if stimulus source
            F = repmat(n2.f, 1, size(DEN, 2));
        else
            F = (abs(NUM1).*n1.f(X1i) + abs(NUM2).*n1.f(X2i) + DEN.*n2.f(Zi)) ...
                ./(abs(NUM1) + abs(NUM2) + DEN);
            con.epsn = n1.e.^((NUM1+NUM2+DEN-2)/2); % used in zdot
        end
        con.epsc = epsilon.^((NUM1+NUM2+DEN-2)/2);  % used in cdot
        con.nType = 4;
        
    case 'active' % Full series of active nonlinearities
        F = (2*F1.*F2)./(F1+F2);
        con.nType = 5;
        
    case 'all2freq' % Full series of resonant monomials
        F = (2*F1.*F2)./(F1+F2);
        con.nType = 6;
        
    case 'allfreq' % Full series of resonant monomials including conjugates
        F = (2*F1.*F2)./(F1+F2);
        con.nType = 7;
        
end

%% Scaling factors

if n2.nFspac==2 % log spacing
    con.F = F;
    con.w = w.*n2.f;
else
    F = ones(size(F));
    con.w = w;
end

%% Scale learning parameters

% if any(kappa) % JCK: replaces if sum(sum(kappa))
%     con.learn = 1;
% else
%     con.learn = 0;
% end

con.learn = learn;

if n2.nFspac==2 % log spacing
    con.lambda = lambda.*F;
    con.mu1 = mu1.*F;
    con.mu2 = mu2.*F;
    con.kappa = kappa.*F;
else
    con.lambda = lambda;
    con.mu1 = mu1;
    con.mu2 = mu2;
    con.kappa = kappa;
end
con.e = epsilon;

%% Connection State and Initial Conditions
%       C0: initial conditions
%        C: connection matrix (instantaneous)
%        t: times saved in 3D memory matrix
%       C3: state memory (3D matrix: frequency x frequency x time)

if isempty(C)
    if ~con.learn
        error('Connection values must be specified as C (3rd input argument) if not learning')
    end
    A0 = zeros(size(F));
    A = spontAmp(real(con.lambda(1,1)), real(con.mu1(1,1)), real(con.mu2(1,1)), con.e);
    A0 = A0+min(A);
    A0 = A0.*(1 +.01*randn(size(A0)));
    theta0 = randn(size(A0));
    con.C0 = A0.*exp(1i*2*pi*theta0);
else
    if isscalar(C)
        con.C0 = C*ones(size(F));
    elseif ~isequal(size(C), size(F))
        error('Incorrect size of C (3rd input argument)')
    else
        con.C0 = C;
    end
end

con.C = con.C0;

%% Masking
%      This process operates on the initial condition and parameter
%      matrices to hange the way things work (i.e.
%
%      mask = ~eye(size(C));          % Don't learn self-connections
%                                       This could be inverted guassian
%                                       too

mask = 1;
if strcmpi(con.type, '2freq') && con.no11
    mask = ~(con.NUM == 1 & con.DEN == 1);
elseif strcmpi(con.type, '3freq')
    mask = con.mask;
elseif con.source == con.target && con.nSourceClass == 2
    mask = ~eye(size(con.C));          % ... Won't learn self-connections
end

con.lambda = con.lambda .* mask;
con.mu1 = con.mu1 .* mask;
con.mu2 = con.mu2 .* mask;
con.kappa = con.kappa .* mask;
con.C0    = con.C0    .* mask;
con.C     = con.C     .* mask;

%% Return network n2
n = n2;
n.con{con.id} = con;

%% Add index to learnList if learn
if con.learn
    n.learnList = [n.learnList con.id];
end

varargin = {@zdot, @cdot, s, n};

    model.zfun = varargin{1};
    model.cfun = varargin{2};
    ind = 3;



%% Get stimuli first to get time info
stimListAll = []; % list of all stimulus id's
stimList = []; % list of stimulus id's, only those used as source
netList = []; % list of network id's
fs = [];
dt = [];
ts = [];
t = [];

for v = ind:length(varargin)
    temp = varargin{v};
    if temp.nClass == 1
        
        sid = temp.id;
        temp.N = size(temp.x, 1);   % in case extra channel was added
        temp.z = temp.x(:,1);   % initialize current state z
        model.s{sid} = temp;
        stimListAll = [stimListAll sid];

            fs = temp.fs;
            dt = temp.dt;
            ts = temp.ts;
            t = temp.t;
    end
end

model.fs           = fs;
model.dt           = dt;
model.t            = t;

%% Get networks and set initial conditions

for v = ind:length(varargin)
    temp = varargin{v};
    if temp.nClass == 2        
        nid = temp.id;
        model.n{nid} = temp;
        netList = [netList nid];
                         
            Nt = ceil(length(t)/temp.sStep);
            model.n{nid}.t = t(1:temp.sStep:length(t));
            model.n{nid}.Z = single(zeros(length(temp.z), Nt));
            model.n{nid}.Z(:,1) = temp.z0;
        
        for cx = temp.learnList
            con = temp.con{cx};
                Nt = ceil(length(t)/con.sStep);
                model.n{nid}.con{cx}.t = t(1:con.sStep:length(t));
                model.n{nid}.con{cx}.C3 = single(zeros(size(con.C,1), size(con.C,2), Nt));
                model.n{nid}.con{cx}.C3(:,:,1) = con.C0;
        end
    end
end

netList = sort(netList);
stimListAll = sort(stimListAll);

%% Check if all connections are valid and at least one network gets stimulus
for j = netList
    for k = 1:length(model.n{j}.con)
        con = model.n{j}.con{k};
        if con.nSourceClass == 1
                stimList = [stimList con.source];            
        end
    end
end

model.stimListAll = stimListAll;
model.stimList = sort(stimList);
model.netList = netList;


%% Cast everything as complex and single

for nx = model.netList
    model.n{nx}.z0 = castCS(model.n{nx}.z0);
    model.n{nx}.z  = castCS(model.n{nx}.z);
    model.n{nx}.Z  = castCS(model.n{nx}.Z);
    model.n{nx}.a  = castCS(model.n{nx}.a);
    model.n{nx}.b1 = castCS(model.n{nx}.b1);
    model.n{nx}.b2 = castCS(model.n{nx}.b2);
    model.n{nx}.e  = castCS(model.n{nx}.e);
    for cx = 1:length(model.n{nx}.con)
        if any(model.n{nx}.learnList) && any(model.n{nx}.learnList == cx)
            model.n{nx}.con{cx}.C0 = castCS(model.n{nx}.con{cx}.C0);
            model.n{nx}.con{cx}.C  = castCS(model.n{nx}.con{cx}.C);
            model.n{nx}.con{cx}.C3 = castCS(model.n{nx}.con{cx}.C3);
        end
        model.n{nx}.con{cx}.w      = castCS(model.n{nx}.con{cx}.w);
        model.n{nx}.con{cx}.lambda = castCS(model.n{nx}.con{cx}.lambda);
        model.n{nx}.con{cx}.mu1    = castCS(model.n{nx}.con{cx}.mu1);
        model.n{nx}.con{cx}.mu2    = castCS(model.n{nx}.con{cx}.mu2);
        model.n{nx}.con{cx}.kappa  = castCS(model.n{nx}.con{cx}.kappa);
        model.n{nx}.con{cx}.e      = castCS(model.n{nx}.con{cx}.e);
    end
end


%% If dotfunc (override) option is empty, then use base dotfunc from network and oscillator-model

if isempty(model.zfun)
    n = varargin{1}; % for now, restrict to only one osc-model
    if strcmp(n.model, 'vdp')
        model.zfun = @zdotv;
    end
    if strcmp(n.model, 'wc')
        model.zfun = @zdotw_sc;
    end
    if strcmp(n.model, 'wce')
        model.zfun = @zdotw_sc;
    end
    if strcmp(n.model, 'hopft')
        model.zfun = @zdotw_sc;
    end
    if strcmp(n.model, 'hopfx')
        model.zfun = @zdotw_sc;
    end
    
end

M = model;

%% Run the network
tic

load('MyColormaps', 'IF_colormap');
circular = IF_colormap;

zfun = M.zfun;
cfun = M.cfun;
ispan = [1 length(M.t)];

%% Error checking
if( ~isa(zfun,'function_handle') )
    error('odeRK4fs: odefun param must be a function handle');
end

step = single(1);
h = single(M.dt);                   % step size
stimList = M.stimList;
netList = M.netList;

%% Display stimulus and initial conditions if dStep > 0
for sx = stimList
    if M.s{sx}.dStep
        stimulusLiveDisplay(M.s{sx}, 0, M.t(1));
    end
end

for nx = netList
    if M.n{nx}.dStep
        networkLiveDisplay(M.n{nx}, 0, M.t(1));
    end
    for cx = M.n{nx}.learnList
        if M.n{nx}.con{cx}.dStep
            connectionLiveDisplay(M.n{nx}.con{cx}, 0, M.t(1), circular);
        end
    end
end

%% Integration loop
for ix = ispan(1) : step : ispan(2)-step
    ind = ix; % time step for which to calculate k1
    
    %% Get Runge-Kutta k-values
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
                ind = ix + step/2; % time step for k2 and k3
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
                ind = ix + step; % time step for k4
            case 4
                for nx = netList
                    M.n{nx}.z = M.n{nx}.zPrev + ...
                        (M.n{nx}.k{1} + 2*M.n{nx}.k{2} + 2*M.n{nx}.k{3} + M.n{nx}.k{4})/6;
                    if M.n{nx}.sStep && ~mod(ix, M.n{nx}.sStep)
                        M.n{nx}.Z(:,ix/M.n{nx}.sStep+1) = M.n{nx}.z;
                    end
                    if M.n{nx}.dStep && ~mod(ix, M.n{nx}.dStep)
                        networkLiveDisplay(M.n{nx}, ix, M.t(ix));
                    end
                    for cx = M.n{nx}.learnList
                        M.n{nx}.con{cx}.C = M.n{nx}.con{cx}.CPrev + ...
                            (M.n{nx}.con{cx}.k{1} + 2*M.n{nx}.con{cx}.k{2} + 2*M.n{nx}.con{cx}.k{3} + M.n{nx}.con{cx}.k{4})/6;
                        if M.n{nx}.con{cx}.sStep && ~mod(ix, M.n{nx}.con{cx}.sStep)
                            M.n{nx}.con{cx}.C3(:,:,ix/M.n{nx}.con{cx}.sStep+1) = M.n{nx}.con{cx}.C;
                        end
                        if M.n{nx}.con{cx}.dStep && ~mod(ix, M.n{nx}.con{cx}.dStep)
                            connectionLiveDisplay(M.n{nx}.con{cx}, ix, M.t(ix), circular);
                        end
                    end
                end
                
                for sx = stimList
                    if M.s{sx}.dStep && ~mod(ix, M.s{sx}.dStep)
                        stimulusLiveDisplay(M.s{sx}, ix, M.t(ix));
                    end
                end
        end
    end
end
toc
