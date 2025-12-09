% Plot ROI histograms as a 3D PDF waterfall for one video.
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
histIdx = 1:nF;
edges = linspace(0,1000,nBins+1);
binCenters = (edges(1:end-1)+edges(2:end))/2;
binWidth = edges(2)-edges(1);

% z positions
if isfield(vid,'zPos') && numel(vid.zPos)>=nF
    zValsAll = vid.zPos(:);
else
    zValsAll = (1:nF).';
end
zSel = zValsAll(histIdx);

H = nan(nBins,nF);
meanInt = nan(nF,1);
for k = 1:nF
    img = frames{k};
    if isempty(img), continue; end
    mask = getRoiMask(vid, size(img));
    pix = double(img(mask));
    counts = histcounts(pix, edges);
    areaCounts = sum(counts)*binWidth;
    if areaCounts>0, counts = counts./areaCounts; end % PDF
    H(:,k) = counts(:);
    meanInt(k) = mean(pix);
end

% reference z*: peak of mean intensity
zRef = zSel(end);
finMean = isfinite(zSel) & isfinite(meanInt);
if nnz(finMean)>=3
    [~, idxMax] = max(meanInt(finMean));
    zRef = parabolicPeak(zSel(finMean), meanInt(finMean), idxMax);
end

fig = figure('Name','ROI histograms 3D (PDF)','Color','w');
set(fig,'Units','normalized','Position',[0 0 1 0.7]);
axw = axes(fig); hold(axw,'on');
baseCmap = autumn(256);
cols = baseCmap(round(linspace(1,size(baseCmap,1), nF)),:);
for k = 1:nF
    [yStep,zStep] = stairs(binCenters, H(:,k));
    xStep = (zSel(k)-zRef)*ones(size(yStep));
    plot3(axw, xStep, yStep, zStep, 'LineWidth', 2, 'Color', cols(k,:));
    xPoly = [xStep; flipud(xStep)];
    yPoly = [yStep; flipud(yStep)];
    zPoly = [zStep; zeros(size(zStep))];
    fill3(axw, xPoly, yPoly, zPoly, cols(k,:), 'FaceAlpha',0.15,'EdgeColor','none');
end
view(axw,45,30);
xlabel(axw,'$z$ ($\mu$m)','Interpreter','latex');
ylabel(axw,'Pixel Intensity','Interpreter','latex');
zlabel(axw,'PDF','Interpreter','latex');
ylim(axw,[edges(1) edges(end)]);
box(axw,'on'); grid(axw,'off');
pbaspect(axw,[10 3 1]);

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

function zPeak = parabolicPeak(zVals, yVals, idxMax)
if nargin<3 || isempty(idxMax), [~,idxMax]=max(yVals); end
zPeak = zVals(idxMax);
if idxMax<=1 || idxMax>=numel(zVals), return; end
z0=zVals(idxMax-1); z1=zVals(idxMax); z2=zVals(idxMax+1);
y0=yVals(idxMax-1); y1=yVals(idxMax); y2=yVals(idxMax+1);
if ~all(isfinite([y0 y1 y2])), return; end
den = (y0 - 2*y1 + y2);
if den==0, return; end
zPeak = z1 + 0.5*((y0 - y2)/den)*(z2 - z0)/2;
end
