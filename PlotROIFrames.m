% Script: plot three representative frames (first, middle, last) with ROI overlay and ROI intensity histogram.
% Edit matFile and vidIndex as needed.

matFile = "\\Actnem\all_protocols_and_methods\XA_Reports_DataAnalysis_Literature\1_RAW_VIDEOS\CONFOCAL\3D_FORMATION_ACTIVE_NEMATICS\20251121_3DACTIVENEMATICS_with_without_ATP\20251121_3DACTIVENEMATICS_ATP\all_videos_roi.mat";
vidIndex = 1; % change if you want another video
nBins = 100;

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

% Plot all three frames in a single aligned figure
fig = figure('Name','ROI frames + histograms','Color','w');
set(fig,'Units','normalized','Position',[0 0 1 1]); % fullscreen
tiledlayout(fig,numel(frameIdx),2,'TileSpacing','compact','Padding','compact');
allPix = {}; % collect ROI pixels to unify histogram axes
for k = 1:numel(frameIdx)
    fi = frameIdx(k);
    img = vid.frames{fi};
    if isempty(img), continue; end
    mask = getRoiMask(vid, size(img));
    roiPixels = double(img(mask));
    allPix{end+1} = roiPixels; %#ok<AGROW>

    ax1 = nexttile((k-1)*2+1); hold(ax1,'on');
    imagesc(ax1, img);
    colormap(ax1, gray);
    axis(ax1,'image'); axis(ax1,'off');
    title(ax1, sprintf('Frame %d', fi));
    % ROI boundary
    B = bwboundaries(mask);
    for b = 1:numel(B)
        plot(ax1, B{b}(:,2), B{b}(:,1), 'r-', 'LineWidth', 1);
    end
    % ROI label
    text(ax1, 5, 15, 'ROI', 'Color','r','FontSize',12,'FontWeight','bold');

    ax2 = nexttile((k-1)*2+2); hold(ax2,'on');
    histogram(ax2, roiPixels, nBins, 'EdgeColor','none','FaceColor',[0.2 0.6 0.9]);
    xlabel(ax2,'Intensity (a.u.)');
    ylabel(ax2,'Count');
    title(ax2,'ROI intensity histogram');
    box(ax2,'on');
end

% unify histogram axes
if ~isempty(allPix)
    allPixCat = cell2mat(allPix(:));
    xmin = min(allPixCat); xmax = max(allPixCat);
    histAxes = findobj(fig,'Type','axes','-not','Tag','legend');
    histAxes = histAxes(arrayfun(@(h) any(h.XLabel.String=="Intensity (a.u.)"), histAxes));
    for h = histAxes'
        xlim(h,[xmin xmax]);
    end
end

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
