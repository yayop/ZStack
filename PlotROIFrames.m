% Script: plot three representative frames (first, middle, last) with ROI overlay and ROI intensity histogram.
% Edit matFile and vidIndex as needed.

matFile = "\\Actnem\all_protocols_and_methods\XA_Reports_DataAnalysis_Literature\1_RAW_VIDEOS\CONFOCAL\3D_FORMATION_ACTIVE_NEMATICS\20251121_3DACTIVENEMATICS_with_without_ATP\20251121_3DACTIVENEMATICS_ATP\all_videos_roi.mat";
vidIndex = 140; % change if you want another video
nBins = 30;

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

% z positions
if isfield(vid,'zPos') && numel(vid.zPos) >= nF
    zValsAll = vid.zPos(:);
else
    zValsAll = (1:nF).';
end
zSel = zValsAll;

% Compute z* as peak of mean ROI intensity
meanIntAll = nan(nF,1);
for k = 1:nF
    imgk = vid.frames{k};
    if isempty(imgk), continue; end
    maskk = getRoiMask(vid, size(imgk));
    meanIntAll(k) = mean(double(imgk(maskk)));
end
zRef = zSel(end);
finMean = isfinite(zSel) & isfinite(meanIntAll);
if nnz(finMean) >= 3
    [~, idxMax] = max(meanIntAll(finMean));
    zRef = parabolicPeak(zSel(finMean), meanIntAll(finMean), idxMax);
elseif nnz(finMean) >= 1
    idxList = find(finMean);
    [~, idxMax] = max(meanIntAll(finMean));
    zRef = zSel(idxList(idxMax));
end

% Output directory: same folder as this script (repo)
scriptDir = fileparts(mfilename('fullpath'));
outDir = scriptDir;
if ~exist(outDir,'dir')
    mkdir(outDir);
end

% Compute global min/max over selected frames for consistent scaling
imgMin = inf; imgMax = -inf;
for k = 1:numel(frameIdx)
    fi = frameIdx(k);
    img = vid.frames{fi};
    if isempty(img), continue; end
    imgMin = min(imgMin, double(min(img(:))));
    imgMax = max(imgMax, double(max(img(:))));
end
imgRange = imgMax - imgMin;
if imgRange <= 0
    imgRange = 1;
end

for k = 1:numel(frameIdx)
    fi = frameIdx(k);
    img = vid.frames{fi};
    if isempty(img), continue; end
    mask = getRoiMask(vid, size(img));
    
    % ROI boundary on full frame (no overlay in export)
    % use fixed scaling based on global min/max from selected frames and boost brightness slightly
    scaleImg = (double(img) - imgMin) / max(imgRange, eps);
    scaleImg = min(max(scaleImg + 0.05, 0), 1); % small brightness boost
    imgUint = im2uint8(scaleImg);
    rgb = repmat(imgUint, [1 1 3]);

    % Save frame with ROI overlay as tif, name includes z shifted by z*
    zShift = zSel(fi) - zRef;
    outName = fullfile(outDir, sprintf('frame_%03d_z_%0.2f.tif', fi, zShift));
    imwrite(rgb, outName);
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

function zPeak = parabolicPeak(zVals, yVals, idxMax)
% Subpixel peak from max point and its neighbors using a parabola fit.
if nargin < 3 || isempty(idxMax)
    [~, idxMax] = max(yVals);
end
zPeak = zVals(idxMax);
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
