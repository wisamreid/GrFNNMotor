disp('initializing workspace');

basePath = '/N/u/madcanda/Mason/grangerCausalityAnalysis/';
dataPath = '/N/u/madcanda/Mason/data/asdf/';
asdfUtilsPath = '/N/u/madcanda/Mason/asdf_utils/ASDF_utils/';
mvgcPath = '/N/u/madcanda/Mason/grangerCausalityAnalysis/mvgc_v1.0/';

addpath(dataPath);
addpath((asdfUtilsPath));
addpath(genpath(mvgcPath));


a = 1;
save('/N/u/madcanda/Mason/grangerCausalityAnalysis/results/test.mat','a');

load('asdf.mat');
N = asdf_raw{end}(1); % number of neurons

% summing up activity at each time bin across neurons 
fs = [];
for i=1:50000:asdf_raw{end}(2)
    cutasdf = ASDFChooseTime(asdf_raw,i,i+50000-1);
    [raster,~] = ASDFToSparse(cutasdf);
    raster = full(raster);
    fs = [fs sum(raster)];
end
clear raster;
clear cutasdf;

% filter
disp('filtering to create mask');
load('spikeBurstMaskingFilter.mat');
fsf = filter(Hd,fs);
% adjust for shift caused by filter
delay = mean(grpdelay(Hd));
fsf = [fsf(delay:end) fsf(1:delay-1)];

% threshold
fsf(fsf>=0.75) = 1;
fsf(fsf<0.75) = 0;

% mask has been prepared
disp('mask has been prepared');
mask = fsf;
save('/N/u/madcanda/Mason/grangerCausalityAnalysis/results/mask.mat','mask');
clear fsf;

% create burst raster
disp('creating burst raster by applying mask to data');
burstRaster = [];
for i=1:N
    dummyAsdf{1,1} = asdf_raw{i};
    dummyAsdf{2,1} = 1;
    dummyAsdf{3,1} = [1,asdf_raw{end}(2)];

    [raster , ~] = ASDFToSparse(dummyAsdf);
    raster = full(raster);
    %apply mask
    raster = raster(find(mask == 1));
    % append to burst raster
    burstRaster = [burstRaster; raster];
end
clear asdf_raw;

% convolve burstRaster to get continous signal
disp('convolving burst raster with gaussian to get continuous signal');
burstSignal = gausConvRaster(burstRaster,1,size(burstRaster,2));
save('/N/u/madcanda/Mason/grangerCausalityAnalysis/results/burstRaster.mat','burstRaster');
clear burstRaster;
save('/N/u/madcanda/Mason/grangerCausalityAnalysis/results/burstSignal.mat','burstSignal');

% Granger Analysis
disp('granger analyzing continuous signal');
%[F,cd,sig] = rasterToMVGC(burstSignal);
X = burstSignal;
clear burstSignal;
disp('\n\nGranger Causality Analysis--');
disp('setting parameters');
ntrials   = 1;     % number of trials

nobs      = size(X,2);   % number of observations per trial

regmode   = 'OLS';  % VAR model estimation regression mode ('OLS', 'LWR' or empty for default)
icregmode = 'LWR';  % information criteria regression mode ('OLS', 'LWR' or empty for default)

morder    = 'AIC';  % model order to use ('actual', 'AIC', 'BIC' or supplied numerical value)
momax     = 20;     % maximum model order for model order estimation

acmaxlags = 1000;   % maximum autocovariance lags (empty for automatic calculation)

tstat     = '';     % statistical test for MVGC:  'F' for Granger's F-test (default) or 'chi2' for Geweke's chi2 test
alpha     = 0.05;   % significance level for significance test
mhtc      = 'FDR';  % multiple hypothesis test correction (see routine 'significance')

fs        = 200;    % sample rate (Hz)
fres      = [];     % frequency resolution (empty for automatic calculation)

seed      = 0;      % random seed (0 for unseeded)


morder = momax;
%% VAR model estimation
% Estimate VAR model of selected order from data.
disp('Estimating VAR model');
ptic('\n*** tsdata_to_var... ');
[A,SIG] = tsdata_to_var(X,morder,regmode);
ptoc;

% Check for failed regression

assert(~isbad(A),'VAR estimation failed');

% NOTE: at this point we have a model and are finished with the data! - all
% subsequent calculations work from the estimated VAR parameters A and SIG.

%% AutoCovariance calculation
% The autocovariance sequence drives many Granger causality calculations (see
% next section). Now we calculate the autocovariance sequence G according to the
% VAR model, to as many lags as it takes to decay to below the numerical
% tolerance level, or to acmaxlags lags if specified (i.e. non-empty).
disp('Computing autocovariance using var model');
ptic('*** var_to_autocov... ');
[G,info] = var_to_autocov(A,SIG,acmaxlags);
ptoc;

% The above routine does a LOT of error checking and issues useful diagnostics.
% If there are problems with your data (e.g. non-stationarity, colinearity,
% etc.) there's a good chance it'll show up at this point - and the diagnostics
% may supply useful information as to what went wrong. It is thus essential to
% report and check for errors here.

var_info(info,true); % report results (and bail out on error)

%% Time domain Granger causality analysis 
% Calculate time-domain pairwise-conditional causalities - this just requires
% the autocovariance sequence.
disp('Time domain Granger Causality Analysis');
ptic('*** autocov_to_pwcgc... ');
F = autocov_to_pwcgc(G);
ptoc;

% Check for failed GC calculation

assert(~isbad(F,false),'GC calculation failed');

% Significance test using theoretical null distribution, adjusting for multiple
% hypotheses.
disp('Significance testing');
pval = mvgc_pval(F,morder,nobs,ntrials,1,1,nvars-2,tstat); % take careful note of arguments!
sig  = significance(pval,alpha,mhtc);

% Plot time-domain causal graph, p-values and significance.
% 
% figure(2); clf;
% subplot(1,3,1);
% plot_pw(F);
% title('Pairwise-conditional GC');
% subplot(1,3,2);
% plot_pw(pval);
% title('p-values');
% subplot(1,3,3);
% plot_pw(sig);
% title(['Significant at p = ' num2str(alpha)])

% For good measure we calculate Seth's causal density (cd) measure - the mean
% pairwise-conditional causality. We don't have a theoretical sampling
% distribution for this.

cd = mean(F(~isnan(F)));


save('/N/u/madcanda/Mason/grangerCausalityAnalysis/results/GCAresults.mat','F');
save('/N/u/madcanda/Mason/grangerCausalityAnalysis/results/causalDensity.mat','cd');
save('/N/u/madcanda/Mason/grangerCausalityAnalysis/results/GCAsig.mat','sig');
