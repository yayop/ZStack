function ZStack_IntensityVsZ(matFile)
%ZSTACK_INTENSITYVSZ Plot mean ROI intensity vs Z for each video.
%   ZStack_IntensityVsZ() loads all_videos_roi.mat (exported by ZStackGUI)
%   from the current folder, or prompts for a file if not found. Each video
%   curve is color-coded by acquisition start time (absTime) if available,
%   otherwise by video index/name.
%
%   ZStack_IntensityVsZ(matFile) loads the specified MAT file.

if nargin < 1 || isempty(matFile)
    defaultFile = fullfile(pwd,'all_videos_roi.mat');
    if exist(defaultFile,'file')
        matFile = defaultFile;
    else
        [f,p] = uigetfile('*.mat','Select all_videos_roi.mat');
        if isequal(f,0), return; end
        matFile = fullfile(p,f);
    end
end

if ~exist(matFile,'file')
    error('File not found: %s', matFile);
end

S = load(matFile);
if ~isfield(S,'roiData')
    error('roiData not found in %s', matFile);
end
roiData = S.roiData;
if ~isstruct(roiData)
    error('roiData is not a struct.');
end

nVids = numel(roiData);
if nVids == 0
    warning('No videos in roiData.');
    return;
end

% Compute color map based on absolute start times if available
[colors, cbLabel] = computeColors(roiData);
figure('Name','Mean ROI intensity vs Z','Color','w');
hold on;
plotted = false(1,nVids);

for v = 1:nVids
    vid = roiData(v);
    if ~isfield(vid,'frames') || isempty(vid.frames)
        continue;
    end
    frames = vid.frames;
    nF = numel(frames);
    if nF == 0
        continue;
    end
    meanVals = nan(nF,1);
    for k = 1:nF
        img = frames{k};
        if isempty(img)
            meanVals(k) = NaN;
        else
            meanVals(k) = mean(double(img(:)),'omitnan');
        end
    end
    zVals = [];
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
    plot(zVals, meanVals, 'Color', colors(v,:), 'LineWidth', 1.5);
    plotted(v) = true;
end

if ~any(plotted)
    warning('No usable frames to plot.');
    return;
end

xlabel('Z position');
ylabel('Mean intensity (ROI)');
title('Mean ROI intensity vs Z');
box on;
colormap(colors);
cb = colorbar;
cb.Label.String = cbLabel;
cb.TickLabelInterpreter = 'none';
if numel(cb.Ticks) ~= size(colors,1)
    % adjust ticks for continuous scale
    cb.Ticks = linspace(cb.Limits(1), cb.Limits(2), min(6, size(colors,1)));
end
if strcmp(cbLabel,'Video index')
    cb.Ticks = 1:nVids;
    cb.TickLabels = arrayfun(@(v)safeName(v), roiData, 'UniformOutput', false);
else
    % show min/max times
    try
        times = extractAbsStart(roiData);
        if ~isempty(times) && any(~isnat(times))
            cb.TickLabels = cellstr(datestr(linspace(min(times(~isnat(times))), max(times(~isnat(times))), numel(cb.Ticks))));
        end
    catch
    end
end

hold off;

end

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

function [colors, label] = computeColors(roiData)
n = numel(roiData);
times = extractAbsStart(roiData);
if isempty(times) || all(isnat(times))
    colors = parula(n);
    label = 'Video index';
    return;
end
label = 'Abs start time';
% map times to [0,1] and use abyss if available
valid = ~isnat(times);
mins = min(times(valid));
maxs = max(times(valid));
if mins == maxs
    tnorm = zeros(n,1);
else
    tnorm = (seconds(times - mins)) ./ seconds(maxs - mins);
end
cmap = getAbyssColormap(256);
idx = 1 + round(tnorm*(size(cmap,1)-1));
colors = cmap(idx,:);
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
            if isdatetime(times) && ~isnat(t0)
                % unify timezones by converting to UTC then removing tz
                try
                    if ~isempty(t0.TimeZone)
                        t0 = datetime(t0,'TimeZone','UTC');
                    end
                    times(i) = datetime(t0,'TimeZone','');
                catch
                    times(i) = t0;
                end
            end
        elseif isnumeric(t) && ~isempty(t)
            % maybe seconds; treat as offset from 0
            times(i) = datetime(t(1), 'ConvertFrom','posixtime');
        end
    end
end
end

function cmap = getAbyssColormap(n)
%GETABYSSCOLORMAP Try to use abyss if installed; fall back to parula.
if nargin < 1, n = 256; end
cmap = [];
if exist('abyss','file') == 2
    try
        cmap = feval('abyss', n);
    catch
        try
            cmap = colormap('abyss');
        catch
        end
    end
else
    try
        cmap = colormap('abyss');
    catch
    end
end
if isempty(cmap)
    cmap = parula(n);
end
% ensure correct length
if size(cmap,1) ~= n
    xo = linspace(0,1,n);
    xi = linspace(0,1,size(cmap,1));
    cmap = interp1(xi, cmap, xo, 'linear');
end
end
