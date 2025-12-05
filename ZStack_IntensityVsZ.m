% Script: plots mean ROI intensity vs Z for each video in all_videos_roi.mat
% Configure the source file here:
matFile = "\\Actnem\all_protocols_and_methods\XA_Reports_DataAnalysis_Literature\1_RAW_VIDEOS\CONFOCAL\3D_FORMATION_ACTIVE_NEMATICS\20251121_3DACTIVENEMATICS_with_without_ATP\20251121_3DACTIVENEMATICS_noATP\all_videos_roi.mat"; % edit if needed
%matFile = "\\Actnem\all_protocols_and_methods\XA_Reports_DataAnalysis_Literature\1_RAW_VIDEOS\CONFOCAL\3D_FORMATION_ACTIVE_NEMATICS\20251121_3DACTIVENEMATICS_with_without_ATP\20251121_3DACTIVENEMATICS_ATP\all_videos_roi.mat";
%matFile = "\\actnem\all_protocols_and_methods\XA_Reports_DataAnalysis_Literature\1_RAW_VIDEOS\CONFOCAL\3D_FORMATION_ACTIVE_NEMATICS\20251129_3DFORMATION_ACTIVENEMATICS\noATP\all_videos_roi.mat";
if ~exist(matFile,'file')
    [f,p] = uigetfile('*.mat','Select all_videos_roi.mat');
    if isequal(f,0), error('File not found.'); end
    matFile = fullfile(p,f);
end

S = load(matFile);
if ~isfield(S,'roiData')
    error('roiData not found in %s', matFile);
end
roiData = S.roiData;
if ~isstruct(roiData) || isempty(roiData)
    error('roiData is empty or not a struct.');
end
%%
nFrames = arrayfun(@(r) numel(r.frames), roiData);
nVids = numel(roiData);
[colors, cbLabel, relTimes, relSpan, baseCmap] = computeColors(roiData);
if nVids == 0
    error('No videos found after loading roiData.');
end
% Subsample curves to at most 14, approximately equispaced
maxCurves = min(14, nVids);
idxPlot = unique(round(linspace(1, nVids, maxCurves)));
roiDataPlot = roiData(idxPlot);
colorsPlot = colors(idxPlot,:);
relTimesPlot = relTimes;
if ~isempty(relTimesPlot)
    relTimesPlot = relTimesPlot(idxPlot);
end
fig = figure('Name','Mean ROI intensity vs Z','Color','w');
set(fig,'Units','normalized','Position',[0 0 1 0.6]);
tiledlayout(fig,1,6,'TileSpacing','compact','Padding','compact');
ax1 = nexttile(1); hold(ax1,'on');

% Reference z0 from latest curve (gaussian peak if possible)
refZ = 0;
% (No z-shift applied; use raw metadata z positions)

minZall = inf; maxZall = -inf;
zMaxList = nan(nVids,1); % now stores centroid z
tList = duration(zeros(nVids,1),0,0); % store as duration
fitA = nan(nVids,1);
fitB = nan(nVids,1);
fitMu = nan(nVids,1);
fitSigma = nan(nVids,1);
zStore = cell(nVids,1);
yStore = cell(nVids,1);
% reference z using the latest time (if available)
if ~isempty(relTimes)
    [~, refVid] = max(relTimes);
    [refMean, refZvals] = computeMeanZ(roiData(refVid));
    if ~isempty(refMean) && ~isempty(refZvals)
        [refZ, ~] = gaussianPeak(refZvals, refMean);
    end
else
    refZ = 0;
end
for v = 1:nVids
    vid = roiData(v);
    [meanVals, zVals] = computeMeanZ(vid);
    if isempty(meanVals), continue; end
    % shift by reference z*
    zVals = zVals - refZ;
    % Ensure monotonic z for integration/peak finding
    [zVals, order] = sort(zVals);
    meanVals = meanVals(order);
    minZall = min(minZall, min(zVals));
    maxZall = max(maxZall, max(zVals));
    % Peak via gaussian fit (fallback to parabolic if needed)
    [zMark, yMark] = gaussianPeak(zVals, meanVals);
    zMaxList(v) = zMark;
    zStore{v} = zVals;
    yStore{v} = meanVals;
    if ~isempty(relTimes)
        tList(v) = relTimes(v);
    end
    % Fit Gaussian model B + A*exp(-(x-mu)^2/(4*sigma^2))
    zFit = zVals(:);
    yFit = meanVals(:);
    nfit = min(numel(zFit), numel(yFit));
    zFit = zFit(1:nfit);
    yFit = yFit(1:nfit);
    fin = isfinite(zFit) & isfinite(yFit);
    nGood = nnz(fin);
    if nGood >= 4
        zFit = zFit(fin);
        yFit = yFit(fin);
        A0 = max(yFit) - min(yFit);
        B0 = min(yFit);
        mu0 = zMark;
        sigma0 = max(std(zFit), eps);
        [~, imaxFit] = max(yFit);
        peakZ = zFit(imaxFit);
        seeds = [
            A0, B0, mu0, sigma0;
            A0, B0, mu0, max(range(zFit)/4, eps);
            A0, B0, peakZ, max(range(zFit)/6, eps)
            ];
        done = false;
        for si = 1:size(seeds,1)
            sp = seeds(si,:);
            opts = fitoptions('Method','NonlinearLeastSquares', ...
                'StartPoint',sp, ...
                'Lower',[0 -Inf min(zFit) 0], ...
                'Upper',[Inf Inf max(zFit) Inf], ...
                'Display','Off', ...
                'Robust','Bisquare');
            ft = fittype('B + A*exp(-((x-mu).^2)/(4*sigma^2))', ...
                'independent','x','coefficients',{'A','B','mu','sigma'}, ...
                'options',opts);
            try
                res = fit(zFit(:), yFit(:), ft, opts);
                fitA(v) = res.A;
                fitB(v) = res.B;
                fitMu(v) = res.mu;
                fitSigma(v) = res.sigma;
                done = true;
                break;
            catch
                % try next seed
            end
        end
        if ~done
            % fallback to simple estimates
            fitA(v) = A0;
            fitB(v) = B0;
            fitMu(v) = mu0;
            fitSigma(v) = sigma0;
        end
    elseif nGood >= 2
        zFit = zFit(fin);
        yFit = yFit(fin);
        fitA(v) = max(yFit) - min(yFit);
        fitB(v) = min(yFit);
        fitMu(v) = zMark;
        fitSigma(v) = max([std(zFit), range(zFit)/4, eps]);
    end
    % Plot only subset
    if ~ismember(v, idxPlot), continue; end
    vidIdxPlot = find(idxPlot==v,1,'first');
    scatter(ax1, zVals, meanVals, 30, 'MarkerFaceColor', colorsPlot(vidIdxPlot,:), 'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 0.9,'LineWidth',0.1);
    % Mark peak point
    % Overlay fitted Gaussian (only for plotted subset)
    if ~isnan(fitA(v)) && ~isnan(fitB(v)) && ~isnan(fitMu(v)) && ~isnan(fitSigma(v)) && fitSigma(v) > 0
        zFine = linspace(min(zVals), max(zVals), 200);
        yFine = fitB(v) + fitA(v) .* exp(-((zFine - fitMu(v)).^2) ./ (4*fitSigma(v).^2));
        plot(ax1, zFine, yFine, 'k--', 'LineWidth', 2);
    end
end

% Minimal legend: one marker (data) and one line (fit)
legendColor = baseCmap(1,:);
demoScatter = scatter(ax1, -inf, -inf, 30, 'MarkerFaceColor',legendColor, ...
    'MarkerEdgeColor', [0 0 0], 'DisplayName','Data','Visible','off');
demoLine = plot(ax1, [-inf -inf], [-inf -inf], 'k--', 'LineWidth', 2, 'DisplayName','Fit','Visible','off');
legend(ax1,[demoScatter,demoLine],{'Data','Fit'},'Location','northwest','Box','on','AutoUpdate','off','FontWeight','bold');

xlabel(ax1,'$z$ ($\mu$m)','Interpreter','latex','FontSize',17);
ylabel(ax1,'$\langle I \rangle$','Interpreter','latex','FontSize',17);
set(ax1,'FontSize',12);
axis(ax1,'square');
pbaspect(ax1,[1 1 1]);
% Fixed x-limits and ticks after shifting by z*
xlim(ax1,[-20 60]);
xticks(ax1,[-20 0 20 40 60]);
box(ax1,'on');
colormap(ax1, baseCmap);
caxis(ax1,[0 relSpan]);
cb = colorbar(ax1);
cb.Limits = [0 relSpan];
cb.TickLabelInterpreter = 'tex'; % keep ticks horizontal/unstyled
cb.Color = [0 0 0];
cb.EdgeColor = 'k';
% Title on top instead of side label
cb.Title.String = '$t$ (min)';
cb.Title.Interpreter = 'latex';
cb.Title.FontSize = 13;
cb.Title.Color = [0 0 0];
% Three integer ticks across the span
tickVals = linspace(0, relSpan, 3);
cb.Ticks = tickVals;
cb.TickLabels = arrayfun(@(t)sprintf('%d', round(t)), tickVals, 'UniformOutput', false);

validFit = ~isnan(fitA) & ~isnan(fitB) & ~isnan(fitMu) & ~isnan(fitSigma);
tminsAll = minutes(tList);
tmins = tminsAll(validFit);
cFit = colors(validFit,:);
xt = [0 20 40 60 80];

axA = nexttile(2); hold(axA,'on');
scatter(axA, tmins, fitA(validFit), 50, cFit, 'filled','MarkerEdgeColor','none');
xlim(axA,[0 80]); xticks(axA, xt); setAdaptiveY(axA, fitA(validFit));
axis(axA,'square'); pbaspect(axA,[1 1 1]); set(axA,'PlotBoxAspectRatio',[1 1 1]);
xlabel(axA,'$t$ (min)','Interpreter','latex'); ylabel(axA,'$A$','Interpreter','latex');

axB = nexttile(3); hold(axB,'on');
scatter(axB, tmins, fitB(validFit), 50, cFit, 'filled','MarkerEdgeColor','none');
xlim(axB,[0 80]); xticks(axB, xt); setAdaptiveY(axB, fitB(validFit));
axis(axB,'square'); pbaspect(axB,[1 1 1]); set(axB,'PlotBoxAspectRatio',[1 1 1]);
xlabel(axB,'$t$ (min)','Interpreter','latex'); ylabel(axB,'$B$','Interpreter','latex');

axMu = nexttile(4); hold(axMu,'on');
scatter(axMu, tmins, fitMu(validFit), 50, cFit, 'filled','MarkerEdgeColor','none');
[absSlopeMu, yLineMu] = addLinFit(axMu, tmins, fitMu(validFit));
xlim(axMu,[0 80]); xticks(axMu, xt); setAdaptiveY(axMu, [fitMu(validFit); yLineMu(:)]);
axis(axMu,'square'); pbaspect(axMu,[1 1 1]); set(axMu,'PlotBoxAspectRatio',[1 1 1]);
title(axMu, sprintf('$|v| = %.2f~(\\mu$m/min)', absSlopeMu),'Interpreter','latex','Color',[0 0 0],'FontSize',14);
xlabel(axMu,'$t$ (min)','Interpreter','latex','FontSize',12); ylabel(axMu,'$\mu$ ($\mu$m)','Interpreter','latex');

axS = nexttile(5); hold(axS,'on');
scatter(axS, tmins, fitSigma(validFit), 50, cFit, 'filled','MarkerEdgeColor','none');
xlim(axS,[0 80]); xticks(axS, xt); setAdaptiveY(axS, fitSigma(validFit));
axis(axS,'square'); pbaspect(axS,[1 1 1]); set(axS,'PlotBoxAspectRatio',[1 1 1]);
xlabel(axS,'$t$ (min)','Interpreter','latex','FontSize',12); ylabel(axS,'$\sigma$ ($\mu$m)','Interpreter','latex');
set([axA axB axMu axS],'FontSize',12,'Box','on');

hold(ax1,'off'); hold(axA,'off'); hold(axB,'off'); hold(axMu,'off'); hold(axS,'off');

% Sixth subplot: collapsed profiles (I-B)/A vs (z-mu)/sigma
axN = nexttile(6); hold(axN,'on');
nGauss = linspace(-4,4,200);
% Model overlay matches fit form: B + A*exp(-(z-mu)^2/(4 sigma^2)) => exp(-0.25 x^2) when x = (z-mu)/sigma
for v = 1:nVids
    if isnan(fitA(v)) || isnan(fitB(v)) || isnan(fitMu(v)) || isnan(fitSigma(v)), continue; end
    if fitA(v) == 0 || fitSigma(v) <= 0, continue; end
    zv = zStore{v}; yv = yStore{v};
    if isempty(zv) || isempty(yv), continue; end
    zstd = (zv - fitMu(v)) ./ fitSigma(v);
    ystd = (yv - fitB(v)) ./ fitA(v);
    scatter(axN, zstd, ystd, 50, 'MarkerFaceColor', colors(v,:), ...
        'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.6, 'LineWidth', 0.1);
end
plot(axN, nGauss, exp(-0.25*nGauss.^2), 'k--', 'LineWidth', 1.5, 'DisplayName','Gaussian overlay');
xlabel(axN,'$(z-\mu)/\sigma$','Interpreter','latex','FontSize',14);
ylabel(axN,'$(I-B)/A$','Interpreter','latex','FontSize',14);
set(axN,'FontSize',12);
axis(axN,'square'); pbaspect(axN,[1 1 1]); box(axN,'on');

% --- Helpers ------------------------------------------------------------
function name = safeName(vid)
if isfield(vid,'stackName') && ~isempty(vid.stackName)
    name = vid.stackName;
elseif isfield(vid,'stackPath') && ~isempty(vid.stackPath)
    [~,nm,ext] = fileparts(vid.stackPath);
    name = [nm ext];
else
    name = sprintf('Video %d', vid.videoIndex);
end
end

function [colors, label, relTimes, relSpan, baseCmap] = computeColors(roiData)
n = numel(roiData);
times = extractAbsStart(roiData);
if isempty(times) || all(isnat(times))
    baseCmap = winter(256);
    colors = baseCmap(round(linspace(1,size(baseCmap,1), n)),:);
    label = 'Video index';
    relTimes = [];
    relSpan = 0;
    return;
end
label = 'Abs start - min (UTC)';
valid = ~isnat(times);
t0 = min(times(valid));
relTimes = times - t0;
mins = min(relTimes(valid));
maxs = max(relTimes(valid));
relSpan = minutes(maxs - mins);
if relSpan <= 0
    relSpan = eps;
end
if mins == maxs
    tnorm = zeros(n,1);
else
    tnorm = seconds(relTimes - mins) ./ seconds(maxs - mins);
end
baseCmap = winter(256);
idx = 1 + round(tnorm*(size(baseCmap,1)-1));
colors = baseCmap(idx,:);
end

function times = extractAbsStart(roiData)
n = numel(roiData);
times = NaT(n,1);
for i = 1:n
    vid = roiData(i);
    if isfield(vid,'absTime') && ~isempty(vid.absTime)
        t = vid.absTime;
        if iscell(t), try t = t{1}; end; end
        if isdatetime(t)
            t0 = t(1);
            try
                if ~isempty(t0.TimeZone)
                    t0 = datetime(t0,'TimeZone','UTC');
                end
                times(i) = datetime(t0,'TimeZone','');
            catch
                times(i) = t0;
            end
        elseif isnumeric(t) && ~isempty(t)
            times(i) = datetime(t(1), 'ConvertFrom','posixtime');
        end
    end
end
end

function [meanVals, zVals] = computeMeanZ(vid)
meanVals = [];
zVals = [];
if ~isfield(vid,'frames') || isempty(vid.frames), return; end
frames = vid.frames;
nF = numel(frames);
if nF == 0, return; end
meanVals = nan(nF,1);
for k = 1:nF
    img = frames{k};
    if isempty(img)
        meanVals(k) = NaN;
    else
        meanVals(k) = median(double(img(:)),'omitnan');
    end
end
if isfield(vid,'zPos')
    zVals = vid.zPos;
end
if isempty(zVals) || numel(zVals) ~= nF || all(isnan(zVals))
    zVals = 1:nF;
else
    zVals = zVals(:).';
    if numel(zVals) ~= nF
        zVals = zVals(1:min(end,nF));
        meanVals = meanVals(1:numel(zVals));
    end
end
end

function cmap = abyssPalette(n)
% Use installed abyss colormap (assumed to exist)
if nargin < 1, n = 256; end
try
    cmap = feval('abyss', n);
catch
    cmap = colormap('abyss');
end
if size(cmap,1) ~= n
    xo = linspace(0,1,n);
    xi = linspace(0,1,size(cmap,1));
    cmap = interp1(xi, cmap, xo, 'linear');
end
end

function [zPeak, yPeak] = interpMax(zVals, yVals)
%INTERPMAX Smooth interpolate (spline if available, then makima, else pchip) and locate peak.
zVals = zVals(:);
yVals = yVals(:);
if numel(zVals) < 3 || numel(yVals) < 3
    [~, idx] = max(yVals);
    zPeak = zVals(idx);
    yPeak = yVals(idx);
    return;
end
zf = linspace(min(zVals), max(zVals), 500);
try
    yf = spline(zVals, yVals, zf);
catch
    try
        yf = interp1(zVals, yVals, zf, 'makima');
    catch
        yf = pchip(zVals, yVals, zf);
    end
end
[yPeak, idx] = max(yf);
zPeak = zf(idx);
end

function zc = centroidZ(zVals, yVals)
zVals = zVals(:);
yVals = yVals(:);
n = min(numel(zVals), numel(yVals));
zVals = zVals(1:n);
yVals = yVals(1:n);
fin = ~isnan(zVals) & ~isnan(yVals);
zVals = zVals(fin);
yVals = yVals(fin);
if isempty(zVals) || all(yVals<=0)
    zc = NaN;
    return;
end
zc = sum(zVals(:).*yVals(:)) ./ sum(yVals(:));
end

function [zPeak, yPeak] = gaussianPeak(zVals, yVals)
zVals = zVals(:); yVals = yVals(:);
n = min(numel(zVals), numel(yVals));
zVals = zVals(1:n); yVals = yVals(1:n);
fin = ~isnan(zVals) & ~isnan(yVals);
zVals = zVals(fin); yVals = yVals(fin);
if numel(zVals) < 3 || all(yVals<=0)
    [yPeak, idx] = max(yVals);
    zPeak = zVals(idx);
    return;
end
[~, imax] = max(yVals);
win = max(1, imax-2):min(numel(zVals), imax+2);
zSeg = zVals(win);
ySeg = yVals(win);
if numel(zSeg) < 3 || all(ySeg<=0)
    [yPeak, idx] = max(yVals);
    zPeak = zVals(idx);
    return;
end
% Initial guesses
a0 = max(ySeg) - min(ySeg);
d0 = min(ySeg);
b0 = zSeg(ySeg==max(ySeg));
c0 = std(zSeg);
gfun = @(p,z) p(1).*exp(-((z-p(2)).^2)./(2*p(3)^2)) + p(4);
p0 = [a0, b0(1), max(c0, eps), d0];
lb = [0, min(zSeg), eps, -Inf];
ub = [Inf, max(zSeg), Inf, Inf];
try
    opts = optimset('Display','off');
    p = lsqcurvefit(gfun, p0, zSeg, ySeg, lb, ub, opts);
    zPeak = p(2);
    yPeak = gfun(p, zPeak);
catch
    % fallback parabola
    p2 = polyfit(zSeg, ySeg, 2);
    if p2(1) ~= 0
        zVert = -p2(2)/(2*p2(1));
        yVert = polyval(p2, zVert);
    else
        [yVert, idx] = max(ySeg);
        zVert = zSeg(idx);
    end
    zPeak = zVert;
    yPeak = yVert;
end
end

function [absSlope, yLine] = addLinFit(ax, tvals, yvals)
absSlope = NaN; yLine = [];
fin = isfinite(tvals) & isfinite(yvals);
if nnz(fin) < 2, return; end
t = tvals(fin); y = yvals(fin);
out = isoutlier(y);
mask = ~out;
if nnz(mask) < 2, mask = fin; end
t = t(mask); y = y(mask);
p = polyfit(t, y, 1);
absSlope = abs(p(1));
tLine = linspace(min(t), max(t), 100);
yLine = polyval(p, tLine);
plot(ax, tLine, yLine, 'k--', 'LineWidth', 3);
end

function setAdaptiveY(ax, ydata)
ydata = ydata(:);
ydata = ydata(isfinite(ydata));
if isempty(ydata), return; end
ymin = min(ydata); ymax = max(ydata);
span = ymax - ymin;
if span == 0, span = max(abs(ymin),1); end
pad = 0.05*span;
ylim(ax, [ymin - pad, ymax + pad]);
axis(ax,'square');
pbaspect(ax,[1 1 1]);
set(ax,'PlotBoxAspectRatio',[1 1 1]);
% keep same x-limits after aspect adjustments
xlim(ax,[0 85]);
end
