%% 
% Script: plot ROI frames (imshow) and a 3D waterfall of ROI histograms across slices.
% Edit matFile and vidIndex as needed.

matFile = "\\Actnem\all_protocols_and_methods\XA_Reports_DataAnalysis_Literature\1_RAW_VIDEOS\CONFOCAL\3D_FORMATION_ACTIVE_NEMATICS\20251121_3DACTIVENEMATICS_with_without_ATP\20251121_3DACTIVENEMATICS_ATP\all_videos_roi.mat";
vidIndex = 140; % change to select a different video
nBins = 200; % more bins for finer histogram
makePerCurveFigures = true; % set true to open a small figure per slice

if ~exist(matFile,'file')
    [f,p] = uigetfile('*.mat','Select all_videos_roi.mat');
    if isequal(f,0), error('File not found.'); end
    matFile = fullfile(p,f);
end

S = load(matFile);
if ~isfield(S,'roiData') || ~isstruct(S.roiData)
    error('roiData not found in %s', matFile);
end
roiData = S.roiData;
if vidIndex < 1 || vidIndex > numel(roiData)
    error('vidIndex out of range. There are %d videos.', numel(roiData));
end
vid = roiData(vidIndex);
if ~isfield(vid,'frames') || isempty(vid.frames)
    error('Selected video has no frames.');
end

nF = numel(vid.frames);
frameIdx = unique([1, round(nF/2), nF]); % first, middle, last
nSel = numel(frameIdx);
histIdx = 1:nF; % use all frames for histogram view
nHist = numel(histIdx);

%% Show selected frames with ROI overlay (imshow, no contrast stretch)
fig1 = figure('Name','Selected frames with ROI','Color','w');
set(fig1,'Units','normalized','Position',[0 0 1 0.5]);
tiledlayout(fig1,1,numel(frameIdx),'TileSpacing','compact','Padding','compact');
% Compute global intensity bounds across selected frames
allImgVals = [];
for k = 1:numel(frameIdx)
    img = vid.frames{frameIdx(k)};
    if isempty(img), continue; end
    allImgVals = [allImgVals; double(img(:))]; %#ok<AGROW>
end
imin = min(allImgVals); imax = max(allImgVals);
for k = 1:numel(frameIdx)
    fi = frameIdx(k);
    img = vid.frames{fi};
    if isempty(img), continue; end
    mask = getRoiMask(vid, size(img));
    nexttile; hold on;
    imshow(img, [imin imax], 'InitialMagnification','fit');
    title(sprintf('Frame %d', fi));
    B = bwboundaries(mask);
    for b = 1:numel(B)
        plot(B{b}(:,2), B{b}(:,1), 'r-', 'LineWidth', 1);
    end
end

%% Plot ROI histograms as 3D lines across slices
edges = linspace(0, 1000, nBins+1); % fixed intensity window
binCenters = (edges(1:end-1)+edges(2:end))/2;
binWidth = edges(2)-edges(1);
H = nan(nBins, nHist);
meanInt = nan(nHist,1);
sigmaInt = nan(nHist,1);
medInt = nan(nHist,1);
meanHist = nan(nHist,1);
sigmaHist = nan(nHist,1);
poiLambdaFit = nan(nHist,1); poiAfit = nan(nHist,1);
% z positions for selected slices
if isfield(vid,'zPos') && numel(vid.zPos) >= nF
    zValsAll = vid.zPos(:);
else
    zValsAll = (1:nF).';
end
% reference z = peak of median ROI intensity across slices
roiMed = nan(nF,1);
for k = 1:nF
    imgk = vid.frames{k};
    if isempty(imgk), continue; end
    maskk = getRoiMask(vid, size(imgk));
    roiVec = double(imgk(maskk));
roiMed(k) = median(roiVec);
end
tmp = roiMed;
tmp(~isfinite(tmp)) = -inf;
[~, idxMaxMed] = max(tmp);
if all(~isfinite(tmp))
    idxMaxMed = 1;
end
zSel = zValsAll(histIdx);
for k = 1:nHist
    idx = histIdx(k);
    img = vid.frames{idx};
    if isempty(img), continue; end
    mask = getRoiMask(vid, size(img));
    roiPixels = double(img(mask));
    counts = histcounts(roiPixels, edges);
    areaCounts = sum(counts)*binWidth;
    if areaCounts > 0
        H(:,k) = counts ./ areaCounts; % normalize to PDF
    else
        H(:,k) = counts;
    end
    meanInt(k) = mean(roiPixels);
    sigmaInt(k) = std(roiPixels);
    medInt(k) = median(roiPixels);
    [meanHist(k), sigmaHist(k)] = weightedStats(binCenters, H(:,k));
end
% reference z*: where I(z) is maximum (parabolic peak around the max)
zRef = zValsAll(end);
finMean = isfinite(zSel) & isfinite(meanInt);
if nnz(finMean) >= 3
    [~, idxMax] = max(meanInt(finMean));
    zRef = parabolicPeak(zSel(finMean), meanInt(finMean), idxMax);
elseif nnz(finMean) >= 1
    idxList = find(finMean);
    [~, idxMax] = max(meanInt(finMean));
    zRef = zSel(idxList(idxMax));
end

fig2 = figure('Name','ROI histograms 3D','Color','w');
set(fig2,'Units','normalized','Position',[0 0 1 0.7]);
axw = axes(fig2); hold(axw,'on');
baseCmap = autumn(256); % match main script
cols = baseCmap(round(linspace(1,size(baseCmap,1), nHist)),:);
for k = 1:nHist
    % stair/step representation so bin plateaus are visible
    [yStep, zStep] = stairs(binCenters, H(:,k));
    xStep = (zSel(k)-zRef)*ones(size(yStep));
    plot3(axw, xStep, yStep, zStep, 'LineWidth', 2, 'Color', cols(k,:));
    % area under curve (transparent)
    xPoly = [xStep; flipud(xStep)];
    yPoly = [yStep; flipud(yStep)];
    zPoly = [zStep; zeros(size(zStep))];
    fill3(axw, xPoly, yPoly, zPoly, cols(k,:), 'FaceAlpha', 0.15, 'EdgeColor', 'none');
    % Poisson fit (LSQ on histogram PDF)
    yData = H(:,k);
        xData = binCenters(:);
        good = isfinite(xData) & isfinite(yData);
        if any(good)
            pSeed = poissonSeed(meanHist(k), max(yData(good)));
            [lamFit, Apoi] = fitPoissonLSQ(xData(good), yData(good), pSeed);
            if ~isnan(lamFit) && lamFit>0
                poiVals = poissonPDFVals(binCenters, lamFit, Apoi);
                plot3(axw, xStep(1)*ones(size(binCenters)), binCenters, poiVals, 'k--', 'LineWidth', 1.5);
                poiLambdaFit(k) = lamFit; poiAfit(k) = Apoi;
            end
        end
    % mark mean intensity on histogram
    [~, idxMeanBin] = min(abs(binCenters - meanInt(k)));
    plot3(axw, xStep(1), binCenters(idxMeanBin), H(idxMeanBin,k), 'o', ...
        'MarkerFaceColor', cols(k,:), 'MarkerEdgeColor', 'k', 'MarkerSize', 10, 'LineWidth', 0.8);
end
view(axw, 45, 30);
set(axw,'XDir','normal'); % standard z-shift axis
xlabel(axw,'$z$ ($\mu$m)','Interpreter','latex');
ylabel(axw,'Pixel Intensity','Interpreter','latex');
zlabel(axw,'PDF','Interpreter','latex');
ylim(axw,[edges(1) edges(end)]);
box(axw,'on'); grid(axw,'off');
pbaspect(axw,[10 3 1]); % equalized box aspect for the combined view

%% Optional: per-slice figure (histogram + fits)
if makePerCurveFigures
    for k = 1:nHist
        figk = figure('Name',sprintf('Hist slice %d', histIdx(k)),'Color','w');
        set(figk,'Units','normalized','Position',[0.2 0.2 0.35 0.3]);
        axk = axes(figk); hold(axk,'on');
        stairs(axk, binCenters, H(:,k), 'Color', cols(k,:), 'LineWidth', 1.5);
        area(axk, binCenters, H(:,k), 'FaceColor', cols(k,:), 'FaceAlpha',0.15, 'EdgeColor','none');
        [~, yMedBin] = min(abs(binCenters - medInt(k)));
        scatter(axk, binCenters(yMedBin), H(yMedBin,k), 60, cols(k,:), 'filled','MarkerEdgeColor','k','LineWidth',0.8);
        xlabel(axk,'Pixel Intensity');
        ylabel(axk,'PDF');
        xlim(axk,[0 1200]);
        % semilogy scaling for individual plots (log only on Y)
        set(axk,'YScale','log');
        yPos = H(:,k);
        yPos = yPos(yPos>0);
        if isempty(yPos)
            yMin = 1; yMax = 1;
        else
            yMin = min(yPos); yMax = max(yPos);
        end
        ylim(axk,[yMin*0.5, yMax*1.5]);
        set(axk,'FontSize',10);
        box(axk,'on');
    end
end

% Overlay of all histograms (PDF)
fig3 = figure('Name','Overlay ROI histograms (PDF)','Color','w');
set(fig3,'Units','normalized','Position',[0 0 1 0.5]);
axOver = axes(fig3); hold(axOver,'on');
for k = 1:nHist
    plot(axOver, binCenters, H(:,k), 'Color', cols(k,:), 'LineWidth', 1.2);
end
set(axOver,'YScale','log');
xlabel(axOver,'Pixel Intensity');
ylabel(axOver,'PDF');
xlim(axOver,[0 1200]);
allPos = H(H>0);
if ~isempty(allPos)
    ylim(axOver,[min(allPos)*0.5, max(allPos)*1.5]);
end
set(axOver,'FontSize',11);
box(axOver,'on');

% --- Helpers ------------------------------------------------------------
function mask = getRoiMask(vid, sz)
% Try to recover a ROI mask; fall back to full frame.
if isfield(vid,'roiMask') && ~isempty(vid.roiMask)
    m = vid.roiMask;
    if iscell(m), m = m{1}; end
    mask = logical(m);
    return;
end
if isfield(vid,'roiPoly') && ~isempty(vid.roiPoly)
    poly = vid.roiPoly;
    if iscell(poly), poly = poly{1}; end
    mask = poly2mask(poly(:,1), poly(:,2), sz(1), sz(2));
    return;
end
% default: full frame
mask = true(sz(1), sz(2));
end

function zPeak = parabolicPeak(zVals, yVals, idxMax)
% Subpixel peak from max point and its neighbors using a parabola fit.
if nargin < 3 || isempty(idxMax)
    [~, idxMax] = max(yVals);
end
zPeak = zVals(idxMax);
% need neighbors
if idxMax <= 1 || idxMax >= numel(zVals)
    return;
end
z0 = zVals(idxMax-1); z1 = zVals(idxMax); z2 = zVals(idxMax+1);
y0 = yVals(idxMax-1); y1 = yVals(idxMax); y2 = yVals(idxMax+1);
if ~all(isfinite([y0 y1 y2]))
    return;
end
den = (y0 - 2*y1 + y2);
if den == 0
    return;
end
zPeak = z1 + 0.5*((y0 - y2)/den)*(z2 - z0)/2;
end

function [m, s] = weightedStats(x, w)
% weighted mean/std (unbiased) for histogram-style data
if nargin < 2 || isempty(w)
    w = ones(size(x));
end
mask = isfinite(x) & isfinite(w) & w>=0;
% Guard sizes (pad or trim weights to match x)
if numel(w) ~= numel(x)
    w = w(:);
    if numel(w) > numel(x)
        w = w(1:numel(x));
    else
        w = [w; zeros(numel(x)-numel(w),1)];
    end
end
x = x(:); w = w(:);
mask = mask(1:numel(x));
x = x(mask); w = w(mask);
if isempty(x) || sum(w)==0
    m = NaN; s = NaN;
    return;
end
m = sum(w.*x)/sum(w);
varw = sum(w.*(x-m).^2)/max(sum(w)-1,1);
s = sqrt(varw);
end

function y = poissonPDFVals(x, lambda, A)
% Poisson PMF evaluated at rounded x and scaled by A
y = zeros(size(x));
pos = x >= 0;
if any(pos)
    xp = round(x(pos));
    y(pos) = A .* exp(-lambda) .* (lambda .^ xp) ./ exp(gammaln(xp+1));
end
end

function seed = poissonSeed(meanX, amp)
if isnan(meanX) || meanX<=0
    seed = [1, max(amp, eps)];
else
    seed = [meanX, max(amp, eps)];
end
end

function [lambdaFit, AFit] = fitPoissonLSQ(x, y, seed)
lambdaFit = NaN; AFit = NaN;
if numel(x) < 3 || numel(y) < 3
    return;
end
if nargin < 3 || numel(seed) < 2
    seed = [max(mean(x), eps), max(y)];
end
logLam0 = log(max(seed(1), eps));
logA0 = log(max(seed(2), eps));
obj = @(p) sum( ( y - poissonPDFVals(x, exp(p(1)), exp(p(2))) ).^2 );
opts = optimset('Display','off');
p0 = [logLam0, logA0];
try
    p = fminsearch(obj, p0, opts);
    lambdaFit = exp(p(1));
    AFit = exp(p(2));
catch
    % leave NaNs
end
end
