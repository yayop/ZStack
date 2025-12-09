% Plot corner plots (posteriors) from saved EMG results.
resultFile = "emg_results_vid140.mat"; % update as needed
sliceList = []; % e.g., [1 10 20]; empty => all slices with samples

if ~exist(resultFile,'file')
    [f,p] = uigetfile('*.mat','Select emg_results_*.mat');
    if isequal(f,0), error('Result file not found'); end
    resultFile = fullfile(p,f);
end
R = load(resultFile);
if ~isfield(R,'emgSamples') || isempty(R.emgSamples)
    error('No emgSamples found in %s', resultFile);
end
samples = R.emgSamples;
if isempty(sliceList)
    sliceList = find(cellfun(@(c) ~isempty(c), samples));
end
labels = {'\mu','\sigma','\tau','A'};
for idx = sliceList(:)'
    samp = samples{idx};
    if isempty(samp) || size(samp,2)~=4
        warning('No samples for slice %d; skipping', idx);
        continue;
    end
    fig = figure('Name',sprintf('Corner slice %d', idx),'Color','w');
    set(fig,'Units','normalized','Position',[0.1 0.1 0.7 0.7]);
    plotCorner(samp, labels);
end

% --- Corner plot helper ---
function plotCorner(samples, labels)
n = size(samples,2);
tiled = tiledlayout(n,n,'TileSpacing','compact','Padding','compact');
for i = 1:n
    for j = 1:n
        ax = nexttile(tiled, (i-1)*n + j);
        hold(ax,'on');
        if i==j
            histogram(ax, samples(:,j), 40, 'FaceColor',[0.2 0.4 0.8],'EdgeColor','none');
            xlabel(ax, labels{j},'Interpreter','latex');
        elseif j < i
            xi = samples(:,j); yi = samples(:,i);
            if exist('ksdensity','file')
                try
                    [f,xiGrid,yiGrid] = ksdensity([xi yi], 'NumPoints',50);
                    fmat = reshape(f, numel(xiGrid), numel(yiGrid));
                    contourf(ax, xiGrid, yiGrid, fmat, 6, 'LineStyle','none');
                    colormap(ax, 'parula');
                catch
                    scatter(ax, xi, yi, 5, 'filled','MarkerFaceAlpha',0.3);
                end
            else
                scatter(ax, xi, yi, 5, 'filled','MarkerFaceAlpha',0.3);
            end
            xlabel(ax, labels{j},'Interpreter','latex');
            ylabel(ax, labels{i},'Interpreter','latex');
        else
            axis(ax,'off');
        end
        set(ax,'FontSize',9);
        box(ax,'on');
    end
end
end
