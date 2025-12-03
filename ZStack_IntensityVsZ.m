% Script: plots mean ROI intensity vs Z for each video in all_videos_roi.mat
% Configure the source file here:
%matFile = "\\Actnem\all_protocols_and_methods\XA_Reports_DataAnalysis_Literature\1_RAW_VIDEOS\CONFOCAL\3D_FORMATION_ACTIVE_NEMATICS\20251121_3DACTIVENEMATICS_with_without_ATP\20251121_3DACTIVENEMATICS_noATP\all_videos_roi.mat"; % edit if needed
matFile = "\\Actnem\all_protocols_and_methods\XA_Reports_DataAnalysis_Literature\1_RAW_VIDEOS\CONFOCAL\3D_FORMATION_ACTIVE_NEMATICS\20251121_3DACTIVENEMATICS_with_without_ATP\20251121_3DACTIVENEMATICS_ATP\all_videos_roi.mat";
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
medianCount = median(nFrames);
keepIdx = (nFrames == medianCount);
roiData = roiData(keepIdx);
nVids = numel(roiData);
if nVids == 0
    error('No videos matched the median frame count filter.');
end
[colors, cbLabel, relTimes, relSpan, baseCmap] = computeColors(roiData);
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
tiledlayout(fig,1,2,'TileSpacing','compact','Padding','compact');
ax1 = nexttile; hold(ax1,'on');

% Reference z0 from latest curve (parabolic peak)
refZ = 0;
if ~isempty(relTimes)
    [~, refVid] = max(relTimes);
    [refMean, refZvals] = computeMeanZ(roiData(refVid));
    if ~isempty(refMean) && ~isempty(refZvals)
        [refZ, ~] = parabolicPeak(refZvals, refMean);
    end
end

minZall = inf; maxZall = -inf;
zMaxList = nan(nVids,1); % now stores centroid z
tList = duration(zeros(nVids,1),0,0); % store as duration
for v = 1:nVids
    vid = roiData(v);
    [meanVals, zVals] = computeMeanZ(vid);
    if isempty(meanVals), continue; end
    zVals = zVals - refZ;
    % Ensure monotonic z for integration/peak finding
    [zVals, order] = sort(zVals);
    meanVals = meanVals(order);
    minZall = min(minZall, min(zVals));
    maxZall = max(maxZall, max(zVals));
    % Peak via parabolic fit
    [zMark, yMark] = parabolicPeak(zVals, meanVals);
    zMaxList(v) = zMark;
    if ~isempty(relTimes)
        tList(v) = relTimes(v);
    end
    % Plot only subset
    if ~ismember(v, idxPlot), continue; end
    vidIdxPlot = find(idxPlot==v,1,'first');
    plot(ax1, zVals, meanVals, 'Color', colorsPlot(vidIdxPlot,:), 'LineWidth', 0.8);
    scatter(ax1, zVals, meanVals, 18, 'MarkerFaceColor', colorsPlot(vidIdxPlot,:), 'MarkerEdgeColor', [0 0 0], 'MarkerFaceAlpha', 0.9);
    % Mark peak point
    yMark = interp1(zVals, meanVals, zMark, 'linear','extrap');
    scatter(ax1, zMark, yMark, 70, 'p', 'MarkerEdgeColor', [0 0 0], ...
        'MarkerFaceColor', 'y', 'LineWidth', 0.8, 'MarkerFaceAlpha', 1);
end
xline(ax1,0,'--','Color',[0.3 0.3 0.3],'LineWidth',1);

xlabel(ax1,'$z$ ($\mu$ m)','Interpreter','latex','FontSize',16);
ylabel(ax1,'$\langle I \rangle$','Interpreter','latex','FontSize',16);
set(ax1,'FontSize',12);
axis(ax1,'square');
pbaspect(ax1,[1 1 1]);
if isfinite(minZall) && isfinite(maxZall)
    xlim(ax1,[minZall maxZall]);
end
box(ax1,'on');
colormap(ax1,baseCmap);
cb = colorbar(ax1);
cb.Label.String = '$t$ (min)';
cb.Label.Interpreter = 'latex';
cb.TickLabelInterpreter = 'none';
cb.EdgeColor = 'none';
caxis(ax1,[0 relSpan]);
tickVals = [0 41 82];
tickVals = tickVals(tickVals <= relSpan);
if isempty(tickVals)
    tickVals = [0 relSpan];
end
cb.Ticks = tickVals;
cb.TickLabels = arrayfun(@(t)sprintf('%.1f', t), tickVals, 'UniformOutput', false);

% Second subplot: z-peak vs time
ax2 = nexttile; hold(ax2,'on');
validTZ = ~isnan(zMaxList) & ~isnan(minutes(tList));
timeVals = minutes(tList(validTZ));
scatter(ax2, timeVals, zMaxList(validTZ), 70, 'p', ...
    'MarkerFaceColor', 'y', 'MarkerEdgeColor', [0 0 0], 'LineWidth', 0.8, 'MarkerFaceAlpha', 1);
xlabel(ax2,'time, t (min)','Interpreter','latex','FontSize',16);
ylabel(ax2,'$z_{\\max}~(\\mu m)$','Interpreter','latex','FontSize',16);
set(ax2,'FontSize',12);
axis(ax2,'square');
pbaspect(ax2,[1 1 1]);
if ~isempty(timeVals)
    xlim(ax2,[min(timeVals) max(timeVals)]);
end
box(ax2,'on');

hold(ax1,'off'); hold(ax2,'off');

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

function [zPeak, yPeak] = parabolicPeak(zVals, yVals)
zVals = zVals(:);
yVals = yVals(:);
n = min(numel(zVals), numel(yVals));
zVals = zVals(1:n); yVals = yVals(1:n);
fin = ~isnan(zVals) & ~isnan(yVals);
zVals = zVals(fin); yVals = yVals(fin);
if numel(zVals) < 3
    [yPeak, idx] = max(yVals);
    zPeak = zVals(idx);
    return;
end
[~, imax] = max(yVals);
idxRange = max(1, imax-1):min(numel(zVals), imax+1);
zSeg = zVals(idxRange);
ySeg = yVals(idxRange);
if numel(zSeg) < 3
    [yPeak, idx] = max(yVals);
    zPeak = zVals(idx);
    return;
end
p = polyfit(zSeg, ySeg, 2);
if p(1) ~= 0
    zVert = -p(2)/(2*p(1));
    yVert = polyval(p, zVert);
else
    zVert = zSeg(2);
    yVert = ySeg(2);
end
if zVert < min(zSeg) || zVert > max(zSeg)
    [yPeak, idx] = max(ySeg);
    zPeak = zSeg(idx);
else
    zPeak = zVert;
    yPeak = yVert;
end
end
