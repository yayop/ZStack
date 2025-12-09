% Plot per-slice ROI PDFs (semilogy) and an overlay figure.
matFile = "\\Actnem\all_protocols_and_methods\XA_Reports_DataAnalysis_Literature\1_RAW_VIDEOS\CONFOCAL\3D_FORMATION_ACTIVE_NEMATICS\20251121_3DACTIVENEMATICS_with_without_ATP\20251121_3DACTIVENEMATICS_ATP\all_videos_roi.mat";
vidIndex = 140; % change as needed
nBins = 200;

if ~exist(matFile,'file')
    [f,p] = uigetfile('*.mat','Select all_videos_roi.mat');
    if isequal(f,0), error('File not found.'); end
    matFile = fullfile(p,f);
end
S = load(matFile);
roiData = S.roiData;
if vidIndex<1 || vidIndex>numel(roiData), error('vidIndex out of range'); end
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
emgSamples = cell(nF,1); % posterior samples after burn
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
    [params, samp] = fitEMGMCMC(binCenters(:), counts(:), [muW, max(sigW, eps), max(sigW, eps), max(counts)]);
    emgParams(k,:) = params;
    emgSamples{k} = samp;
end

% Individual plots
for k = 1:nF
    % Linear scale plot
    figLin = figure('Name',sprintf('Hist slice %d (lin)', k),'Color','w');
    set(figLin,'Units','normalized','Position',[0.2 0.2 0.35 0.25]);
    axLin = axes(figLin); hold(axLin,'on');
    stairs(axLin, binCenters, H(:,k), 'Color', cols(k,:), 'LineWidth', 1.5);
    area(axLin, binCenters, H(:,k), 'FaceColor', cols(k,:), 'FaceAlpha',0.15, 'EdgeColor','none');
    if all(isfinite(emgParams(k,:))) && emgParams(k,2)>0 && emgParams(k,3)>0
        emgVals = emgPDF(binCenters, emgParams(k,1), emgParams(k,2), emgParams(k,3), emgParams(k,4));
        plot(axLin, binCenters, emgVals, 'k--','LineWidth',1.2);
        gaussComp = emgParams(k,4) * (1./(emgParams(k,2)*sqrt(2*pi))) .* exp(-0.5*((binCenters - emgParams(k,1))./emgParams(k,2)).^2);
        expComp = emgParams(k,4) * (1./emgParams(k,3)) .* exp(-(binCenters - emgParams(k,1))./emgParams(k,3));
        expComp(binCenters<emgParams(k,1)) = 0;
        plot(axLin, binCenters, gaussComp, 'Color',[0.2 0.2 0.2], 'LineStyle',':','LineWidth',1);
        plot(axLin, binCenters, expComp, 'Color',[0 0.5 0], 'LineStyle','-.','LineWidth',1);
    end
    [~, yMedBin] = min(abs(binCenters - medInt(k)));
    scatter(axLin, binCenters(yMedBin), H(yMedBin,k), 60, 'filled','MarkerEdgeColor','k','LineWidth',0.8);
    xlabel(axLin,'Pixel Intensity'); ylabel(axLin,'PDF');
    xlim(axLin,[0 1200]);
    yPos = H(:,k); yPos = yPos(yPos>0);
    if isempty(yPos), yMin=1; yMax=1; else, yMin=min(yPos); yMax=max(yPos); end
    ylim(axLin,[0, max(yPos)*1.1]);
    set(axLin,'FontSize',10); box(axLin,'on');

    % Semilogy plot
    figLog = figure('Name',sprintf('Hist slice %d (log)', k),'Color','w');
    set(figLog,'Units','normalized','Position',[0.2 0.55 0.35 0.25]);
    axLog = axes(figLog); hold(axLog,'on');
    stairs(axLog, binCenters, H(:,k), 'Color', cols(k,:), 'LineWidth', 1.5);
    area(axLog, binCenters, H(:,k), 'FaceColor', cols(k,:), 'FaceAlpha',0.15, 'EdgeColor','none');
    if all(isfinite(emgParams(k,:))) && emgParams(k,2)>0 && emgParams(k,3)>0
        emgVals = emgPDF(binCenters, emgParams(k,1), emgParams(k,2), emgParams(k,3), emgParams(k,4));
        plot(axLog, binCenters, emgVals, 'k--','LineWidth',1.2);
        gaussComp = emgParams(k,4) * (1./(emgParams(k,2)*sqrt(2*pi))) .* exp(-0.5*((binCenters - emgParams(k,1))./emgParams(k,2)).^2);
        expComp = emgParams(k,4) * (1./emgParams(k,3)) .* exp(-(binCenters - emgParams(k,1))./emgParams(k,3));
        expComp(binCenters<emgParams(k,1)) = 0;
        plot(axLog, binCenters, gaussComp, 'Color',[0.2 0.2 0.2], 'LineStyle',':','LineWidth',1);
        plot(axLog, binCenters, expComp, 'Color',[0 0.5 0], 'LineStyle','-.','LineWidth',1);
    end
    scatter(axLog, binCenters(yMedBin), H(yMedBin,k), 60, 'filled','MarkerEdgeColor','k','LineWidth',0.8);
    xlabel(axLog,'Pixel Intensity'); ylabel(axLog,'PDF');
    xlim(axLog,[0 1200]);
    if isempty(yPos), yMin=1; yMax=1; else, yMin=min(yPos); yMax=max(yPos); end
    set(axLog,'YScale','log');
    ylim(axLog,[yMin*0.5, yMax*1.5]);
    set(axLog,'FontSize',10); box(axLog,'on');

    % Posterior plots for parameters
    if ~isempty(emgSamples{k})
        figp = figure('Name',sprintf('Posteriors EMG slice %d',k),'Color','w');
        set(figp,'Units','normalized','Position',[0.25 0.25 0.5 0.5]);
        axp = tiledlayout(figp,2,2,'TileSpacing','compact','Padding','compact');
        labels = {'\mu','\sigma','\tau','A'};
        samp = emgSamples{k};
        for pi = 1:4
            axh = nexttile(axp);
            histogram(axh, samp(:,pi), 40, 'FaceColor',[0.2 0.4 0.8],'EdgeColor','none');
            xlabel(axh, labels{pi},'Interpreter','latex');
            ylabel(axh,'Count');
            set(axh,'FontSize',10);
            box(axh,'on');
        end
    end
end

% Overlay
figO = figure('Name','Overlay ROI PDFs','Color','w');
set(figO,'Units','normalized','Position',[0 0 1 0.5]);
axO = axes(figO); hold(axO,'on');
cols = autumn(max(nF,2)); % match main script palette
for k = 1:nF
    plot(axO, binCenters, H(:,k), 'Color', cols(k,:), 'LineWidth', 1.2);
end
set(axO,'YScale','log');
xlabel(axO,'Pixel Intensity'); ylabel(axO,'PDF');
xlim(axO,[0 1200]);
allPos = H(H>0);
if ~isempty(allPos)
    ylim(axO,[min(allPos)*0.5, max(allPos)*1.5]);
end
set(axO,'FontSize',11); box(axO,'on');

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
% Exponential-modified Gaussian PDF scaled by A
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

function [params, samples] = fitEMGMCMC(x, y, seed)
% Simple Metropolis MCMC to fit EMG (mu,sigma,tau,A); returns MAP params and samples
params = [NaN NaN NaN NaN];
samples = [];
if numel(x) < 3 || numel(y) < 3, return; end
if nargin < 3 || numel(seed) < 4
    seed = [mean(x), std(x), std(x), max(y)];
end
% use only points with y>0
mask = isfinite(x) & isfinite(y) & y>0;
x = x(mask); y = y(mask);
if numel(x) < 3, return; end
logy = log(y);

% initial params in log-space for positives
mu0 = seed(1);
logSig0 = log(max(seed(2), eps));
logTau0 = log(max(seed(3), eps));
logA0 = log(max(seed(4), eps));
p = [mu0, logSig0, logTau0, logA0];

logPost = @(pvec) logLikelihoodEMG(x, logy, pvec);

nIter = 500000; burn = 50000;
step = [0.1, 0.1, 0.1, 0.1]; % proposal std
bestP = p; bestLL = -inf;
samples = zeros(nIter,4);
currLL = logPost(p);
for i = 1:nIter
    pProp = p + step.*randn(1,4);
    llProp = logPost(pProp);
    alpha = exp(llProp - currLL);
    if rand < alpha
        p = pProp;
        currLL = llProp;
    end
    samples(i,:) = p;
    if currLL > bestLL
        bestLL = currLL;
        bestP = p;
    end
end
% take MAP (best) after burn
useSamp = samples(burn+1:end,:);
if isempty(useSamp)
    params = [bestP(1), exp(bestP(2)), exp(bestP(3)), exp(bestP(4))];
    return;
end
llSamp = -inf(size(useSamp,1),1);
for ii = 1:size(useSamp,1)
    llSamp(ii) = logPost(useSamp(ii,:));
end
[~,idx] = max(llSamp);
pMAP = useSamp(idx,:);

params = [pMAP(1), exp(pMAP(2)), exp(pMAP(3)), exp(pMAP(4))];
samples = useSamp;
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
