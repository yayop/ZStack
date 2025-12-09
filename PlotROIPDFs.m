% Plot per-slice ROI PDFs (semilogy) and an overlay figure.
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
edges = linspace(0,1000,nBins+1);
binCenters = (edges(1:end-1)+edges(2:end))/2;
binWidth = edges(2)-edges(1);

if isfield(vid,'zPos') && numel(vid.zPos)>=nF
    zValsAll = vid.zPos(:);
else
    zValsAll = (1:nF).';
end

H = nan(nBins,nF);
medInt = nan(nF,1);
for k = 1:nF
    img = frames{k};
    if isempty(img), continue; end
    mask = getRoiMask(vid, size(img));
    pix = double(img(mask));
    counts = histcounts(pix, edges);
    areaCounts = sum(counts)*binWidth;
    if areaCounts>0, counts = counts./areaCounts; end % PDF
    H(:,k) = counts(:);
    medInt(k) = median(pix);
end

% Individual plots
for k = 1:nF
    figk = figure('Name',sprintf('Hist slice %d', k),'Color','w');
    set(figk,'Units','normalized','Position',[0.2 0.2 0.35 0.3]);
    axk = axes(figk); hold(axk,'on');
    stairs(axk, binCenters, H(:,k), 'LineWidth', 1.5);
    area(axk, binCenters, H(:,k), 'FaceAlpha',0.15, 'EdgeColor','none');
    [~, yMedBin] = min(abs(binCenters - medInt(k)));
    scatter(axk, binCenters(yMedBin), H(yMedBin,k), 60, 'filled','MarkerEdgeColor','k','LineWidth',0.8);
    xlabel(axk,'Pixel Intensity');
    ylabel(axk,'PDF');
    xlim(axk,[0 1200]);
    set(axk,'YScale','log');
    yPos = H(:,k);
    yPos = yPos(yPos>0);
    if isempty(yPos), yMin=1; yMax=1; else, yMin=min(yPos); yMax=max(yPos); end
    ylim(axk,[yMin*0.5, yMax*1.5]);
    set(axk,'FontSize',10); box(axk,'on');
end

% Overlay
figO = figure('Name','Overlay ROI PDFs','Color','w');
set(figO,'Units','normalized','Position',[0 0 1 0.5]);
axO = axes(figO); hold(axO,'on');
cols = autumn(max(nF,2));
for k = 1:nF
    plot(axO, binCenters, H(:,k), 'Color', cols(k,:), 'LineWidth', 1.2);
end
set(axO,'YScale','log');
xlabel(axO,'Pixel Intensity'); ylabel(axO,'PDF');
xlim(axO,[0 1200]);
allPos = H(H>0);
if ~isempty(allPos)
    ylim(axO,[min(allPos)*0.5, max(allPos)*1.5]);
end
set(axO,'FontSize',11); box(axO,'on');

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
