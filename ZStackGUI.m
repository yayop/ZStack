function ZStackGUI()
%ZSTACKGUI Basic ND2 viewer with two axes and right-side controls.
%   ZStackGUI() launches a UI that loads ND2 files via Bio-Formats
%   (place bfmatlab inside the ZStackGUI folder). Left side shows a 1x2
%   layout of axes; the current frame is drawn on the left axis. The right
%   panel hosts load buttons and a selector panel (dropdown + frame slider).

% ---- State -------------------------------------------------------------
app = struct();
app.videos = struct('name',{},'path',{},'stack',{},'frameCount',{},'meta',{});
app.currVideoIdx = 0;
app.currFrame = 1;
app.defaultRoot = '\\Actnem\all_protocols_and_methods\XA_Reports_DataAnalysis_Literature\1_RAW_VIDEOS\CONFOCAL\3D_FORMATION_ACTIVE_NEMATICS';
app.roiHandle = [];
app.roiListeners = event.listener.empty;
app.histBins = 64;

% ---- UI ---------------------------------------------------------------
app.fig = uifigure('Name','ZStack GUI','Position',[100 100 1280 760]);
mainGrid = uigridlayout(app.fig,[1 2]);
mainGrid.ColumnWidth = {'3x','1x'};
mainGrid.RowHeight = {'1x'};

% Left: two axes laid out as subplot(1,2)
leftPanel = uipanel(mainGrid,'Title','View','FontWeight','bold');
leftPanel.Layout.Row = 1;
leftPanel.Layout.Column = 1;
leftGrid = uigridlayout(leftPanel,[1 2]);
leftGrid.ColumnWidth = {'1x','1x'};
leftGrid.RowHeight = {'1x'};
app.axCurrent = uiaxes(leftGrid,'Box','on');
app.axCurrent.Layout.Row = 1; app.axCurrent.Layout.Column = 1;
title(app.axCurrent,'Z Slice');
app.axSecondary = uiaxes(leftGrid,'Box','on');
app.axSecondary.Layout.Row = 1; app.axSecondary.Layout.Column = 2;
title(app.axSecondary,'Histogram');
axis(app.axCurrent,'ij');
app.axSecondary.YDir = 'normal';
axis(app.axSecondary,'equal');
colormap(app.axCurrent,'gray');
colormap(app.axSecondary,'gray');

% Right: controls
rightPanel = uipanel(mainGrid,'Title','Controls','FontWeight','bold');
rightPanel.Layout.Row = 1;
rightPanel.Layout.Column = 2;
rightGrid = uigridlayout(rightPanel,[7 1]);
rightGrid.RowHeight = {30,30,30,30,190,'1x',220};
rightGrid.Padding = [10 10 10 10];
rightGrid.RowSpacing = 8;

app.statusLabel = uilabel(rightGrid,'Text','Ready. Load an ND2 file.','FontWeight','bold');
app.statusLabel.Layout.Row = 1;

app.loadOneBtn = uibutton(rightGrid,'Text','Load ND2','ButtonPushedFcn',@(src,evt)onLoadSingle());
app.loadOneBtn.Layout.Row = 2;

app.loadManyBtn = uibutton(rightGrid,'Text','Load folder','ButtonPushedFcn',@(src,evt)onLoadBatch());
app.loadManyBtn.Layout.Row = 3;

app.infoLabel = uilabel(rightGrid,'Text','No video loaded','HorizontalAlignment','left');
app.infoLabel.Layout.Row = 4;

% ROI panel
roiPanel = uipanel(rightGrid,'Title','Define ROI','FontWeight','bold');
roiPanel.Layout.Row = 5;
roiGrid = uigridlayout(roiPanel,[7 2]);
roiGrid.RowHeight = {28,24,24,24,24,28,24};
roiGrid.ColumnWidth = {70,'1x'};
roiGrid.Padding = [6 6 6 6];
roiGrid.RowSpacing = 4;
roiGrid.ColumnSpacing = 6;

app.drawRoiBtn = uibutton(roiGrid,'Text','Draw ROI','ButtonPushedFcn',@(src,evt)onDrawRoi());
app.drawRoiBtn.Layout.Row = 1;
app.drawRoiBtn.Layout.Column = [1 2];

lblX = uilabel(roiGrid,'Text','X:','HorizontalAlignment','right');
lblX.Layout.Row = 2; lblX.Layout.Column = 1;
lblY = uilabel(roiGrid,'Text','Y:','HorizontalAlignment','right');
lblY.Layout.Row = 3; lblY.Layout.Column = 1;
lblW = uilabel(roiGrid,'Text','Width:','HorizontalAlignment','right');
lblW.Layout.Row = 4; lblW.Layout.Column = 1;
lblH = uilabel(roiGrid,'Text','Height:','HorizontalAlignment','right');
lblH.Layout.Row = 5; lblH.Layout.Column = 1;

app.roiXField = uieditfield(roiGrid,'numeric','Value',0,'Limits',[0 Inf]);
app.roiXField.Layout.Row = 2; app.roiXField.Layout.Column = 2;
app.roiYField = uieditfield(roiGrid,'numeric','Value',0,'Limits',[0 Inf]);
app.roiYField.Layout.Row = 3; app.roiYField.Layout.Column = 2;
app.roiWField = uieditfield(roiGrid,'numeric','Value',0,'Limits',[0 Inf]);
app.roiWField.Layout.Row = 4; app.roiWField.Layout.Column = 2;
app.roiHField = uieditfield(roiGrid,'numeric','Value',0,'Limits',[0 Inf]);
app.roiHField.Layout.Row = 5; app.roiHField.Layout.Column = 2;

app.applyRoiBtn = uibutton(roiGrid,'Text','Apply values','ButtonPushedFcn',@(src,evt)onApplyRoiFromFields());
app.applyRoiBtn.Layout.Row = 6;
app.applyRoiBtn.Layout.Column = [1 2];

lblBins = uilabel(roiGrid,'Text','Bins:','HorizontalAlignment','right');
lblBins.Layout.Row = 7; lblBins.Layout.Column = 1;
app.binField = uieditfield(roiGrid,'numeric','Value',app.histBins,'Limits',[1 Inf],...
    'RoundFractionalValues','on','ValueChangedFcn',@(src,evt)onBinsChanged());
app.binField.Layout.Row = 7; app.binField.Layout.Column = 2;

selectorPanel = uipanel(rightGrid,'Title','Selector','FontWeight','bold');
selectorPanel.Layout.Row = 7;
selectorGrid = uigridlayout(selectorPanel,[5 1]);
selectorGrid.RowHeight = {28,28,22,50,30};
selectorGrid.RowSpacing = 6;
selectorGrid.Padding = [6 6 6 6];

app.videoDropdown = uidropdown(selectorGrid,...
    'Items',{'(no videos loaded)'},...
    'ItemsData',[],...
    'ValueChangedFcn',@(src,evt)onVideoSelected());
app.videoDropdown.Layout.Row = 1;

btnRow = uigridlayout(selectorGrid,[1 2]);
btnRow.Layout.Row = 2;
btnRow.ColumnWidth = {'1x','1x'};
btnRow.RowHeight = {24};
btnRow.ColumnSpacing = 6;
btnRow.Padding = [0 0 0 0];

app.prevVideoBtn = uibutton(btnRow,'Text','Prev video','ButtonPushedFcn',@(src,evt)onPrevVideo());
app.prevVideoBtn.Layout.Row = 1; app.prevVideoBtn.Layout.Column = 1;
app.nextVideoBtn = uibutton(btnRow,'Text','Next video','ButtonPushedFcn',@(src,evt)onNextVideo());
app.nextVideoBtn.Layout.Row = 1; app.nextVideoBtn.Layout.Column = 2;

app.frameLabel = uilabel(selectorGrid,'Text','Frame 0/0','HorizontalAlignment','left');
app.frameLabel.Layout.Row = 3;

app.frameSlider = uislider(selectorGrid,...
    'Limits',[1 2],...
    'Value',1,...
    'MajorTicks',[1 2],...
    'ValueChangingFcn',@(src,evt)onFrameChanged(round(evt.Value)),...
    'ValueChangedFcn',@(src,evt)onFrameChanged(round(src.Value)));
app.frameSlider.Layout.Row = 4;

app.saveBtn = uibutton(selectorGrid,'Text','Save','ButtonPushedFcn',@(src,evt)onSaveCropped());
app.saveBtn.Layout.Row = 5;

% Initial placeholders
resetAxes();

% ---- Nested callbacks -------------------------------------------------
    function onLoadSingle()
        ensureBioFormatsPath();
        startDir = getStartDir();
        [file, path] = uigetfile({'*.nd2','ND2 files (*.nd2)'}, 'Select ND2 file', startDir);
        if isequal(file,0)
            return;
        end
        addVideo(fullfile(path,file));
    end

    function onLoadBatch()
        ensureBioFormatsPath();
        startDir = getStartDir();
        folder = uigetdir(startDir,'Select folder containing ND2 files');
        if isequal(folder,0)
            return;
        end
        files = dir(fullfile(folder,'**','*.nd2'));
        if isempty(files)
            showAlert('No .nd2 files were found in that folder.');
            return;
        end
        filePaths = cellfun(@(f,p)fullfile(f,p), {files.folder}, {files.name}, 'UniformOutput', false);
        [newVideos, nFailed] = loadVideosParallel(filePaths);
        if isempty(newVideos)
            showAlert('No videos were loaded.');
            return;
        end
        integrateVideos(newVideos);
        if nFailed > 0
            app.statusLabel.Text = sprintf('Loaded %d video(s); %d failed.', numel(newVideos), nFailed);
        else
            app.statusLabel.Text = sprintf('Loaded %d video(s).', numel(newVideos));
        end
    end

    function addVideo(filePath)
        [stack, frameCount, meta] = readNd2(filePath);
        if isempty(stack)
            return;
        end
        vid.name = buildDisplayName(filePath);
        vid.path = filePath;
        vid.stack = stack;
        vid.frameCount = frameCount;
        vid.meta = meta;
        app.videos(end+1) = vid; %#ok<AGROW>
        app.currVideoIdx = numel(app.videos);
        app.currFrame = 1;
        clearRoiHandle();
        resetRoiFields();
        updateVideoDropdown();
        updateSlider();
        redrawFrame();
        app.statusLabel.Text = sprintf('Loaded: %s', vid.name);
    end

    function [stack, frameCount, meta] = readNd2(filePath)
        ensureBioFormatsPath();
        stack = [];
        frameCount = 0;
        meta = struct('deltaT',[],'zPos',[],'acqDate','');
        try
            data = bfopen(filePath);
        catch ME
            showAlert(sprintf('Could not open %s: %s', filePath, ME.message));
            return;
        end
        series = data{1};
        frameCount = size(series,1);
        if frameCount == 0
            showAlert('ND2 appears empty.');
            return;
        end
        sample = series{1,1};
        sz = size(sample);
        stack = zeros(sz(1), sz(2), frameCount, 'like', sample);
        for idx = 1:frameCount
            stack(:,:,idx) = series{idx,1};
        end
        meta = readMetadata(filePath, frameCount);
    end

    function onFrameChanged(requestedFrame)
        if isempty(app.videos)
            app.frameSlider.Value = 1;
            return;
        end
        nFrames = app.videos(app.currVideoIdx).frameCount;
        newFrame = min(max(1, requestedFrame), nFrames);
        app.currFrame = newFrame;
        app.frameSlider.Value = newFrame;
        redrawFrame();
    end

    function onVideoSelected()
        if isempty(app.videos) || isempty(app.videoDropdown.ItemsData)
            return;
        end
        app.currVideoIdx = app.videoDropdown.Value;
        app.currFrame = 1;
        clearRoiHandle();
        resetRoiFields();
        updateSlider();
        redrawFrame();
    end

    function onPrevVideo()
        if isempty(app.videos) || isempty(app.videoDropdown.ItemsData)
            return;
        end
        newIdx = max(1, app.currVideoIdx - 1);
        if newIdx == app.currVideoIdx
            return;
        end
        app.videoDropdown.Value = newIdx;
        onVideoSelected();
    end

    function onSaveCropped()
        if isempty(app.videos)
            showAlert('Load a video before saving.');
            return;
        end
        if isempty(app.roiHandle) || ~isvalid(app.roiHandle)
            showAlert('Define an ROI before saving.');
            return;
        end
        rect = round(app.roiHandle.Position);
        rect(3) = max(rect(3),1);
        rect(4) = max(rect(4),1);

        nVids = numel(app.videos);
        roiData = repmat(struct('stackName','','stackPath','','roiPosition',[],'histBins',app.histBins,...
            'frames',[],'deltaT',[],'zPos',[],'acqDate','', 'frameCount',0,'videoIndex',0), nVids,1);
        for v = 1:nVids
            vid = app.videos(v);
            stack = vid.stack;
            sz = size(stack);
            if numel(sz) < 3
                sz(3) = 1;
            end
            w = sz(2); h = sz(1);
            rectAdj = rect;
            rectAdj(1) = min(max(1, rect(1)), w);
            rectAdj(2) = min(max(1, rect(2)), h);
            rectAdj(3) = min(rect(3), w - rectAdj(1) + 1);
            rectAdj(4) = min(rect(4), h - rectAdj(2) + 1);
            rectAdj(3) = max(rectAdj(3),1);
            rectAdj(4) = max(rectAdj(4),1);
            frames = cell(vid.frameCount,1);
            for k = 1:vid.frameCount
                frames{k} = imcrop(stack(:,:,k), rectAdj);
            end
            meta = safeMeta(vid);
            roiData(v).stackName = vid.name;
            roiData(v).stackPath = vid.path;
            roiData(v).roiPosition = rectAdj;
            roiData(v).histBins = app.histBins;
            roiData(v).frames = frames;
            roiData(v).deltaT = meta.deltaT;
            roiData(v).zPos = meta.zPos;
            roiData(v).acqDate = meta.acqDate;
            roiData(v).frameCount = vid.frameCount;
            roiData(v).videoIndex = v;
        end

        [file,path] = uiputfile('*.mat','Save cropped ROI frames', fullfile(getStartDir(), 'all_videos_roi.mat'));
        if isequal(file,0)
            return;
        end
        save(fullfile(path,file),'roiData','-v7.3');
        app.statusLabel.Text = sprintf('Saved %d video(s) ROI frames to %s', nVids, file);
    end

    function onNextVideo()
        if isempty(app.videos) || isempty(app.videoDropdown.ItemsData)
            return;
        end
        newIdx = min(numel(app.videos), app.currVideoIdx + 1);
        if newIdx == app.currVideoIdx
            return;
        end
        app.videoDropdown.Value = newIdx;
        onVideoSelected();
    end

    function onDrawRoi()
        if isempty(app.videos)
            showAlert('Load a video before defining ROI.');
            return;
        end
        clearRoiHandle();
        try
            h = drawrectangle(app.axCurrent,'Color',[0 0.8 1],'FaceAlpha',0.1);
        catch
            showAlert('Could not start ROI drawing.');
            return;
        end
        app.roiHandle = h;
        attachRoiListeners();
        syncRoiFieldsFromHandle();
    end

    function onApplyRoiFromFields()
        if isempty(app.videos)
            showAlert('Load a video before defining ROI.');
            return;
        end
        pos = [app.roiXField.Value, app.roiYField.Value, app.roiWField.Value, app.roiHField.Value];
        pos(3) = max(pos(3), eps);
        pos(4) = max(pos(4), eps);
        if isempty(app.roiHandle) || ~isvalid(app.roiHandle)
            try
                h = drawrectangle(app.axCurrent,'Color',[0 0.8 1],'FaceAlpha',0.1,'Position',pos);
            catch
                showAlert('Could not create ROI with those values.');
                return;
            end
            app.roiHandle = h;
        else
            app.roiHandle.Position = pos;
        end
        attachRoiListeners();
        syncRoiFieldsFromHandle();
        drawHistogram(app.videos(app.currVideoIdx).stack(:,:,app.currFrame));
    end

    function syncRoiFieldsFromHandle()
        if isempty(app.roiHandle) || ~isvalid(app.roiHandle)
            return;
        end
        pos = app.roiHandle.Position;
        app.roiXField.Value = pos(1);
        app.roiYField.Value = pos(2);
        app.roiWField.Value = pos(3);
        app.roiHField.Value = pos(4);
        drawHistogram(app.videos(app.currVideoIdx).stack(:,:,app.currFrame));
    end

    function detachRoiListeners()
        if ~isempty(app.roiListeners)
            for k = 1:numel(app.roiListeners)
                if isvalid(app.roiListeners(k))
                    delete(app.roiListeners(k));
                end
            end
        end
        app.roiListeners = event.listener.empty;
    end

    function attachRoiListeners()
        if isempty(app.roiHandle) || ~isvalid(app.roiHandle)
            app.roiHandle = [];
            detachRoiListeners();
            return;
        end
        detachRoiListeners();
        app.roiListeners(1) = addlistener(app.roiHandle,'MovingROI',@(src,evt)syncRoiFieldsFromHandle());
        app.roiListeners(2) = addlistener(app.roiHandle,'ROIMoved',@(src,evt)syncRoiFieldsFromHandle());
    end

    function clearRoiHandle()
        detachRoiListeners();
        if ~isempty(app.roiHandle) && isvalid(app.roiHandle)
            delete(app.roiHandle);
        end
        app.roiHandle = [];
    end

    function resetRoiFields()
        app.roiXField.Value = 0;
        app.roiYField.Value = 0;
        app.roiWField.Value = 0;
        app.roiHField.Value = 0;
    end

    function onBinsChanged()
        val = max(1, round(app.binField.Value));
        app.histBins = val;
        app.binField.Value = val;
        if ~isempty(app.videos)
            vid = app.videos(app.currVideoIdx);
            drawHistogram(vid.stack(:,:,app.currFrame));
        end
    end

    function redrawFrame()
        if isempty(app.videos)
            resetAxes();
            return;
        end
        vid = app.videos(app.currVideoIdx);
        img = vid.stack(:,:,app.currFrame);
        delete(findall(app.axCurrent,'Type','image'));
        imagesc(app.axCurrent, img);
        axis(app.axCurrent,'image');
        app.axCurrent.XTick = [];
        app.axCurrent.YTick = [];
        title(app.axCurrent, sprintf('Frame %d/%d', app.currFrame, vid.frameCount));
        drawHistogram(img);
        app.frameLabel.Text = sprintf('Frame %d/%d', app.currFrame, vid.frameCount);
        app.infoLabel.Text = sprintf('Video: %s', vid.name);
    end

    function meta = readMetadata(filePath, frameCount)
        meta = struct('deltaT',nan(frameCount,1),'zPos',nan(frameCount,1),'acqDate','');
        try
            r = bfGetReader(filePath);
            r.setSeries(0);
            ms = r.getMetadataStore();
            try
                dtStr = ms.getImageAcquisitionDate(0);
                if ~isempty(dtStr)
                    meta.acqDate = char(dtStr.toString());
                end
            catch
            end
            for idx = 0:frameCount-1
                try
                    t = ms.getPlaneDeltaT(0, idx);
                    if ~isempty(t)
                        meta.deltaT(idx+1) = double(t.value());
                    end
                catch
                end
                try
                    z = ms.getPlanePositionZ(0, idx);
                    if ~isempty(z)
                        meta.zPos(idx+1) = double(z.value());
                    end
                catch
                end
            end
            r.close();
        catch
            % keep defaults
        end
    end

    function resetAxes()
        clearRoiHandle();
        resetRoiFields();
        cla(app.axCurrent); cla(app.axSecondary);
        text(app.axCurrent,0.5,0.5,'Load an ND2','HorizontalAlignment','center','Units','normalized','Color',[0.4 0.4 0.4]);
        axis(app.axCurrent,'off');
        text(app.axSecondary,0.5,0.5,'Define ROI to view histogram','HorizontalAlignment','center','Units','normalized','Color',[0.4 0.4 0.4]);
        app.axSecondary.YDir = 'normal';
        axis(app.axSecondary,'square');
        axis(app.axSecondary,'off');
    end

    function drawHistogram(img)
        cla(app.axSecondary);
        title(app.axSecondary,'Histogram');
        app.axSecondary.YDir = 'normal';
        nb = max(1, round(app.binField.Value));
        app.histBins = nb;
        app.binField.Value = nb;
        if isempty(app.roiHandle) || ~isvalid(app.roiHandle)
            text(app.axSecondary,0.5,0.5,'Define ROI to view histogram','HorizontalAlignment','center','Units','normalized','Color',[0.4 0.4 0.4]);
            axis(app.axSecondary,'square');
            axis(app.axSecondary,'off');
            return;
        end
        try
            mask = createMask(app.roiHandle);
        catch
            mask = [];
        end
        if isempty(mask) || ~any(mask(:))
            text(app.axSecondary,0.5,0.5,'ROI is empty','HorizontalAlignment','center','Units','normalized','Color',[0.4 0.4 0.4]);
            axis(app.axSecondary,'square');
            axis(app.axSecondary,'off');
            return;
        end
        vals = double(img(mask));
        if isempty(vals)
            text(app.axSecondary,0.5,0.5,'ROI is empty','HorizontalAlignment','center','Units','normalized','Color',[0.4 0.4 0.4]);
            axis(app.axSecondary,'square');
            axis(app.axSecondary,'off');
            return;
        end
        if numel(unique(vals)) == 1
            c = numel(vals);
            bar(app.axSecondary, vals(1), c, 'FaceColor',[0.3 0.3 0.9],'EdgeColor','none');
            xlim(app.axSecondary, [vals(1)-0.5, vals(1)+0.5]);
        else
            edges = linspace(min(vals), max(vals), nb+1);
            counts = histcounts(vals, edges);
            centers = edges(1:end-1);
            bar(app.axSecondary, centers, counts, 'FaceColor',[0.3 0.3 0.9],'EdgeColor','none','BarWidth',1);
            xlim(app.axSecondary,[edges(1) edges(end)]);
        end
        axis(app.axSecondary,'tight');
        axis(app.axSecondary,'square');
        box(app.axSecondary,'on');
        app.axSecondary.XTickMode = 'auto';
        app.axSecondary.YTickMode = 'auto';
    end

    function updateSlider()
        if isempty(app.videos)
            app.frameSlider.Limits = [1 2];
            app.frameSlider.MajorTicks = [1 2];
            app.frameSlider.Value = 1;
            app.frameLabel.Text = 'Frame 0/0';
            return;
        end
        n = app.videos(app.currVideoIdx).frameCount;
        app.frameSlider.Limits = [1 max(2,n)];
        ticks = unique(round(linspace(1,n,min(6,n))));
        app.frameSlider.MajorTicks = ticks;
        app.frameSlider.Value = min(max(1, app.currFrame), n);
        app.frameLabel.Text = sprintf('Frame %d/%d', app.currFrame, n);
    end

    function updateVideoDropdown()
        names = {app.videos.name};
        app.videoDropdown.Items = names;
        app.videoDropdown.ItemsData = num2cell(1:numel(names));
        app.videoDropdown.Value = app.currVideoIdx;
    end

    function meta = safeMeta(vid)
        n = vid.frameCount;
        meta = struct('deltaT',nan(n,1),'zPos',nan(n,1),'acqDate','');
        if isfield(vid,'meta') && ~isempty(vid.meta)
            if isfield(vid.meta,'deltaT') && numel(vid.meta.deltaT)>=n
                meta.deltaT = vid.meta.deltaT(:);
            end
            if isfield(vid.meta,'zPos') && numel(vid.meta.zPos)>=n
                meta.zPos = vid.meta.zPos(:);
            end
            if isfield(vid.meta,'acqDate')
                meta.acqDate = vid.meta.acqDate;
            end
        end
    end

    function integrateVideos(newVideos)
        if isempty(newVideos)
            return;
        end
        app.videos = [app.videos, newVideos];
        app.currVideoIdx = numel(app.videos);
        app.currFrame = 1;
        clearRoiHandle();
        resetRoiFields();
        updateVideoDropdown();
        updateSlider();
        redrawFrame();
    end

    function [videosOut, nFailed] = loadVideosParallel(filePaths)
        n = numel(filePaths);
        videosOut = struct('name',{},'path',{},'stack',{},'frameCount',{},'meta',{});
        nFailed = 0;
        if n == 0
            return;
        end
        ensureBioFormatsPath(); % main session path
        usePar = license('test','Distrib_Computing_Toolbox') && ~isempty(ver('parallel'));
        stacks = cell(1,n);
        frameCounts = zeros(1,n);
        names = cell(1,n);
        metas = cell(1,n);
        successes = false(1,n);
        if usePar
            try
                parfor k = 1:n
                    fp = filePaths{k};
                    stack = [];
                    frameCount = 0;
                    nm = '';
                    metaLocal = struct('deltaT',[],'zPos',[],'acqDate','');
                    ok = false;
                    baseDir = fileparts(mfilename('fullpath'));
                    bfDir = fullfile(baseDir,'bfmatlab');
                    if exist(bfDir,'dir')
                        addpath(genpath(bfDir));
                        if exist('bfCheckJavaPath','file')
                            try
                                bfCheckJavaPath(false);
                            catch
                            end
                        end
                    end
                    try
                        data = bfopen(fp);
                        series = data{1};
                        frameCount = size(series,1);
                        if frameCount > 0
                            sample = series{1,1};
                            sz = size(sample);
                            tmp = zeros(sz(1), sz(2), frameCount, 'like', sample);
                            for idx = 1:frameCount
                                tmp(:,:,idx) = series{idx,1};
                            end
                            stack = tmp;
                            [~, b, e] = fileparts(fp);
                            nm = [b e];
                            % metadata
                            metaLocal = struct('deltaT',nan(frameCount,1),'zPos',nan(frameCount,1),'acqDate','');
                            try
                                r = bfGetReader(fp);
                                r.setSeries(0);
                                ms = r.getMetadataStore();
                                try
                                    dtStr = ms.getImageAcquisitionDate(0);
                                    if ~isempty(dtStr)
                                        metaLocal.acqDate = char(dtStr.toString());
                                    end
                                catch
                                end
                                for idx = 0:frameCount-1
                                    try
                                        t = ms.getPlaneDeltaT(0, idx);
                                        if ~isempty(t)
                                            metaLocal.deltaT(idx+1) = double(t.value());
                                        end
                                    catch
                                    end
                                    try
                                        z = ms.getPlanePositionZ(0, idx);
                                        if ~isempty(z)
                                            metaLocal.zPos(idx+1) = double(z.value());
                                        end
                                    catch
                                    end
                                end
                                r.close();
                            catch
                            end
                            ok = true;
                        end
                    catch
                    end
                    stacks{k} = stack;
                    frameCounts(k) = frameCount;
                    names{k} = nm;
                    metas{k} = metaLocal;
                    successes(k) = ok;
                end
            catch
                usePar = false;
            end
        end
        if ~usePar
            for k = 1:n
                fp = filePaths{k};
                stack = [];
                frameCount = 0;
                nm = '';
                metaLocal = struct('deltaT',[],'zPos',[],'acqDate','');
                ok = false;
                try
                    data = bfopen(fp);
                    series = data{1};
                    frameCount = size(series,1);
                    if frameCount > 0
                        sample = series{1,1};
                        sz = size(sample);
                        tmp = zeros(sz(1), sz(2), frameCount, 'like', sample);
                        for idx = 1:frameCount
                            tmp(:,:,idx) = series{idx,1};
                        end
                        stack = tmp;
                        [~, b, e] = fileparts(fp);
                        nm = [b e];
                        metaLocal = readMetadata(fp, frameCount);
                        ok = true;
                    end
                catch
                end
                stacks{k} = stack;
                frameCounts(k) = frameCount;
                names{k} = nm;
                metas{k} = metaLocal;
                successes(k) = ok;
            end
        end
        for k = 1:n
            if successes(k)
                v.name = names{k};
                v.path = filePaths{k};
                v.stack = stacks{k};
                v.frameCount = frameCounts(k);
                v.meta = metas{k};
                videosOut(end+1) = v; %#ok<AGROW>
            end
        end
        nFailed = n - numel(videosOut);
    end

    function ensureBioFormatsPath()
        baseDir = fileparts(mfilename('fullpath'));
        bfDir = fullfile(baseDir,'bfmatlab');
        if exist(bfDir,'dir')
            addpath(genpath(bfDir));
            if exist('bfCheckJavaPath','file')
                try
                    bfCheckJavaPath(false);
                catch
                    % ignore; bfopen will throw if Java path is missing
                end
            end
        end
    end

    function showAlert(msg)
        try
            uialert(app.fig,msg,'Notice');
        catch
            warning(msg);
        end
    end

    function startDir = getStartDir()
        if isfield(app,'defaultRoot') && exist(app.defaultRoot,'dir')
            startDir = app.defaultRoot;
        else
            startDir = pwd;
        end
    end

    function name = buildDisplayName(filePath)
        [~, base, ext] = fileparts(filePath);
        name = [base ext];
    end
end
