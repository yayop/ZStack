function ZStack_IntensityVsZ(matFile)
%ZSTACK_INTENSITYVSZ Plot mean ROI intensity vs Z for each video.
%   ZStack_IntensityVsZ() loads all_videos_roi.mat (exported by ZStackGUI)
%   from the current folder, or prompts for a file if not found. Each video
%   curve is color-coded; a colorbar maps color to video index/name.
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

colors = parula(nVids);
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
cb.Label.String = 'Video index';
cb.Ticks = 1:nVids;
cb.TickLabels = arrayfun(@(v)safeName(v), roiData, 'UniformOutput', false);
cb.TickLabelInterpreter = 'none';

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
