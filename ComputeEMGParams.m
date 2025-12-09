% Compute EMG fits for all slices of a selected video and save results.
% Uses MCMC (Metropolis) and stores MAP params plus thinned samples.

matFile = "\\Actnem\all_protocols_and_methods\XA_Reports_DataAnalysis_Literature\1_RAW_VIDEOS\CONFOCAL\3D_FORMATION_ACTIVE_NEMATICS\20251121_3DACTIVENEMATICS_with_without_ATP\20251121_3DACTIVENEMATICS_ATP\all_videos_roi.mat";
vidIndex = 140;      % video to process
nBins    = 200;      % histogram bins
nIter    = 1000000;  % MCMC iterations
burn     = 100000;   % burn-in
thin     = 1;        % keep every sample after burn (no thinning)
maxStore = 50000;    % cap stored samples per slice

if ~exist(matFile,'file')
    [f,p] = uigetfile('*.mat','Select all_videos_roi.mat');
    if isequal(f,0), error('File not found.'); end
    matFile = fullfile(p,f);
end
S = load(matFile);
roiData = S.roiData;
if vidIndex < 1 || vidIndex > numel(roiData)
    error('vidIndex out of range');
end
vid = roiData(vidIndex);

frames = vid.frames;
nF = numel(frames);
edges = linspace(0,1000,nBins+1);
binCenters = (edges(1:end-1)+edges(2:end))/2;
binWidth = edges(2)-edges(1);

if isfield(vid,'zPos') && numel(vid.zPos)>=nF
    zValsAll = vid.zPos(:);
else
    zValsAll = (1:nF).';
end

H = nan(nBins,nF);
medInt = nan(nF,1);
emgParams = nan(nF,4); % [mu sigma tau A]
emgSamples = cell(nF,1);

parfor k = 1:nF
    img = frames{k};
    if isempty(img), continue; end
    mask = getRoiMask(vid, size(img));
    pix = double(img(mask));
    counts = histcounts(pix, edges);
    areaCounts = sum(counts)*binWidth;
    if areaCounts>0, counts = counts./areaCounts; end % PDF
    H(:,k) = counts(:);
    medInt(k) = median(pix);
    [muW, sigW] = weightedStats(binCenters, counts);
    [params, samp] = fitEMGMCMC(binCenters(:), counts(:), ...
        [muW, max(sigW, eps), max(sigW, eps), max(counts)], ...
        nIter, burn, thin, maxStore);
    emgParams(k,:) = params;
    emgSamples{k} = samp;
end

save(sprintf('emg_results_vid%03d.mat', vidIndex), ...
    'emgParams','emgSamples','H','binCenters','edges','zValsAll','vidIndex','matFile','medInt','nBins','nIter','burn','thin','maxStore','-v7.3');
fprintf('Saved EMG results for vid %d\n', vidIndex);

% --- helpers ---
function mask = getRoiMask(vid, sz)
if isfield(vid,'roiMask') && ~isempty(vid.roiMask)
    m = vid.roiMask; if iscell(m), m = m{1}; end
    mask = logical(m); return;
end
if isfield(vid,'roiPoly') && ~isempty(vid.roiPoly)
    poly = vid.roiPoly; if iscell(poly), poly = poly{1}; end
    mask = poly2mask(poly(:,1), poly(:,2), sz(1), sz(2)); return;
end
mask = true(sz(1), sz(2));
end

function [m, s] = weightedStats(x, w)
if nargin < 2 || isempty(w), w = ones(size(x)); end
mask = isfinite(x) & isfinite(w) & w>=0;
x = x(mask); w = w(mask);
if isempty(x) || sum(w)==0
    m = NaN; s = NaN; return;
end
m = sum(w.*x)/sum(w);
varw = sum(w.*(x-m).^2)/max(sum(w)-1,1);
s = sqrt(varw);
end

function y = emgPDF(x, mu, sigma, tau, A)
sigma = max(sigma, eps);
tau = max(tau, eps);
lambda = 1./tau;
y = zeros(size(x));
pos = isfinite(x);
if any(pos)
    xp = x(pos);
    arg = (sigma.^2 .* lambda) - (xp - mu);
    y(pos) = (A .* lambda./2) .* exp((lambda/2).*(2*mu + lambda.*sigma.^2 - 2*xp)) ...
        .* erfc(arg ./ (sqrt(2).*sigma));
end
end

function [params, samplesOut] = fitEMGMCMC(x, y, seed, nIter, burn, thin, maxStore)
params = [NaN NaN NaN NaN];
samplesOut = [];
if numel(x) < 3 || numel(y) < 3, return; end
if nargin < 3 || numel(seed) < 4
    seed = [mean(x), std(x), std(x), max(y)];
end
if nargin < 4, nIter = 50000; end
if nargin < 5, burn = 5000; end
if nargin < 6, thin = 50; end
if nargin < 7, maxStore = 5000; end

mask = isfinite(x) & isfinite(y) & y>0;
x = x(mask); y = y(mask);
if numel(x) < 3, return; end
logy = log(y);

% initial params (mu, log sigma, log tau, log A)
mu0 = seed(1);
p = [mu0, log(max(seed(2), eps)), log(max(seed(3), eps)), log(max(seed(4), eps))];

logPost = @(pvec) logLikelihoodEMG(x, logy, pvec);

step = [0.1, 0.1, 0.1, 0.1];
bestP = p; bestLL = -inf;
currLL = logPost(p);
store = [];
for i = 1:nIter
    pProp = p + step.*randn(1,4);
    llProp = logPost(pProp);
    alpha = exp(llProp - currLL);
    if rand < alpha
        p = pProp; currLL = llProp;
    end
    if i > burn && mod(i, thin)==0
        if size(store,1) < maxStore
            store = [store; p]; %#ok<AGROW>
        end
    end
    if currLL > bestLL
        bestLL = currLL; bestP = p;
    end
end
if isempty(store)
    store = bestP;
end
params = [bestP(1), exp(bestP(2)), exp(bestP(3)), exp(bestP(4))];
samplesOut = store;
end

function ll = logLikelihoodEMG(x, logy, pvec)
mu = pvec(1); sig = exp(pvec(2)); tau = exp(pvec(3)); A = exp(pvec(4));
if sig<=0 || tau<=0 || A<=0
    ll = -inf; return;
end
yModel = emgPDF(x, mu, sig, tau, A);
if any(~isfinite(yModel)) || any(yModel<=0)
    ll = -inf; return;
end
ll = -sum( (logy - log(yModel)).^2 ); % pseudo log-likelihood in log-space
end
