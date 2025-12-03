% Script: plots mean ROI intensity vs Z for each video in all_videos_roi.mat
% Configure the source file here:
matFile = "\\Actnem\all_protocols_and_methods\XA_Reports_DataAnalysis_Literature\1_RAW_VIDEOS\CONFOCAL\3D_FORMATION_ACTIVE_NEMATICS\20251121_3DACTIVENEMATICS_with_without_ATP\20251121_3DACTIVENEMATICS_noATP\all_videos_roi.mat"; % edit if needed

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
nVids = numel(roiData);
[colors, cbLabel, relTimes] = computeColors(roiData);
figure('Name','Mean ROI intensity vs Z','Color','w');
hold on;

% Reference z0 from latest curve
refZ = 0;
if ~isempty(relTimes)
    [~, refVid] = max(relTimes);
    [refMean, refZvals] = computeMeanZ(roiData(refVid));
    if ~isempty(refMean) && ~isempty(refZvals)
        [~, idxMax] = max(refMean);
        refZ = refZvals(idxMax);
    end
end

for v = 1:nVids
    vid = roiData(v);
    [meanVals, zVals] = computeMeanZ(vid);
    if isempty(meanVals), continue; end
    zVals = zVals - refZ;
    plot(zVals, meanVals, 'Color', colors(v,:), 'LineWidth', 0.8);
    scatter(zVals, meanVals, 18, 'MarkerFaceColor', colors(v,:), 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.9);
end

xlabel('$z~(\\mu m)$','Interpreter','latex','FontSize',16);
ylabel('$\\langle I \\rangle$','Interpreter','latex','FontSize',16);
set(gca,'FontSize',12);
axis square;
pbaspect([1 1 1]);
box on;
colormap(colors);
cb = colorbar;
cb.Label.String = cbLabel;
cb.TickLabelInterpreter = 'none';
if strcmp(cbLabel,'Video index')
    cb.Ticks = 1:nVids;
    cb.TickLabels = arrayfun(@(v)safeName(v), roiData, 'UniformOutput', false);
else
    cb.Label.String = 'time, t (min)';
    relMin = min(relTimes);
    relMax = max(relTimes);
    cb.Ticks = linspace(relMin, relMax, min(6,nVids));
    cb.TickLabels = arrayfun(@(t) sprintf('%.2f', t), cb.Ticks, 'UniformOutput', false);
end

hold off;

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

function [colors, label, relTimes] = computeColors(roiData)
n = numel(roiData);
times = extractAbsStart(roiData);
if isempty(times) || all(isnat(times))
    colors = abyssPalette(n);
    label = 'Video index';
    relTimes = [];
    return;
end
label = 'Abs start - min (UTC)';
valid = ~isnat(times);
t0 = min(times(valid));
relTimes = times - t0;
mins = min(relTimes(valid));
maxs = max(relTimes(valid));
if mins == maxs
    tnorm = zeros(n,1);
else
    tnorm = seconds(relTimes - mins) ./ seconds(maxs - mins);
end
colors = abyssPalette(256);
idx = 1 + round(tnorm*(size(colors,1)-1));
colors = colors(idx,:);
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
        meanVals(k) = mean(double(img(:)),'omitnan');
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
