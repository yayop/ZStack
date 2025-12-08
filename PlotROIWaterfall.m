% Script: plot ROI frames (imshow) and a 3D waterfall of ROI histograms across slices.
% Edit matFile and vidIndex as needed.

matFile = "\\actnem\all_protocols_and_methods\XA_Reports_DataAnalysis_Literature\1_RAW_VIDEOS\CONFOCAL\3D_FORMATION_ACTIVE_NEMATICS\20251129_3DFORMATION_ACTIVENEMATICS\noATP\all_videos_roi.mat";
vidIndex = 1; % change to select a different video
nBins = 50;

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

%% Show selected frames with ROI overlay (imshow, no contrast stretch)
fig1 = figure('Name','Selected frames with ROI','Color','w');
set(fig1,'Units','normalized','Position',[0 0 1 0.5]);
tiledlayout(fig1,1,numel(frameIdx),'TileSpacing','compact','Padding','compact');
for k = 1:numel(frameIdx)
    fi = frameIdx(k);
    img = vid.frames{fi};
    if isempty(img), continue; end
    mask = getRoiMask(vid, size(img));
    nexttile; hold on;
    imshow(img, [], 'InitialMagnification','fit');
    title(sprintf('Frame %d', fi));
    B = bwboundaries(mask);
    for b = 1:numel(B)
        plot(B{b}(:,2), B{b}(:,1), 'r-', 'LineWidth', 1);
    end
end

%% Build waterfall of ROI histograms across slices
edges = linspace(0, double(max(cellfun(@(f) max(f(:)), vid.frames))), nBins+1);
binCenters = (edges(1:end-1)+edges(2:end))/2;
H = nan(nBins, nF);
for k = 1:nF
    img = vid.frames{k};
    if isempty(img), continue; end
    mask = getRoiMask(vid, size(img));
    roiPixels = double(img(mask));
    H(:,k) = histcounts(roiPixels, edges);
end

fig2 = figure('Name','ROI histogram waterfall','Color','w');
set(fig2,'Units','normalized','Position',[0 0 1 0.7]);
% Use bar3-like surface
[Y,Z] = meshgrid(1:nF, binCenters);
surf(Y, Z, H, 'EdgeColor','interp', 'FaceColor','none');
view(45,30);
xlabel('Slice index');
ylabel('Intensity (a.u.)');
zlabel('Count');
colormap(parula);
box on; grid on;

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
