%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Spatial Configuration Analysis of Fly Interactions
%
% Related to Cell manuscript "Male-Male Interactions Shape Mate Selection
% in Drosophila" by Tom Hindmarsh Sten, Rufei Li, Florian Hollunder, Shade
% Eleazer, and Vanessa Ruta
%
% Description:
% This script analyzes spatial configurations and behaviors (singing and 
% wing flicking) in triadic fly interactions. The key objectives include:
%   - Computing and visualizing the spatial relationships between flies
%     during interactions in low-dimensional space
%   - Identifying regions of high behavioral activity (singing, flicking).
%   - Quantifying the likelihood of behaviors occurring in specific spatial
%     configurations.
%
% Key Features:
%   - Calculates pairwise distances and angles between flies in the triad.
%   - Uses t-SNE for dimensionality reduction to identify patterns in 
%     spatial data.
%   - Constructs Gaussian-smoothed heatmaps for behavioral likelihood.
%   - Reconstructs fly positions in a coordinate frame centered on the 
%     female and visualizes spatial configurations of males relative to her.
%   - Accurately identifies high-probability regions for singing and flicking
%     using histogram binning and image processing techniques.
%
% Inputs:
%   - Spatial and behavioral data, with pre-classified singing/flicking
%   events. The loaded variable, allMMF, is a 1x5 structure with 16 fields
%   detailing different behavioral parameters. Flicking and Singing are 
%   classified by JAABA models. For each parameter, there is a 3xNFrames
%   double that tracks each of the three animals. The two males are always
%   rows 1 & 2, while the female is always row 3. 
%
% Outputs:
%   - Heatmaps showing behavioral likelihoods in spatial configurations.
%   - 2D position of animals during which specific behaviors are more
%   likely to occur
%
% Key Dependencies:
%   - Image Processing Toolbox
%   - Circular Statistics Toolbox
%   - Customized function: computeRelativeAngle.m
%
% Binning resolution, Gaussian filter size, and t-SNE parameters can be 
% adjusted for specific datasets.
%
% Author: Tom Hindmarsh Sten
% Date: 12/2/2024
% Affiliation: The Rockefeller University
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
% Load preprocessed dataset, set seed
load('exMMF_CantonS_behaviorForTSNE.mat')
%%
rng(1995)

% Initialize variables
minSongFrames = 30; % Minimum frames to qualify as a song event
minFlickFrames = 3; % Minimum frames to qualify as a flick event
flyCtr = 0; % Counter for unique flies

% Create structures to store results
store.rawPredictors = [];
store.filteredPredictors = [];
store.song = [];
store.flick = [];

% Loop over all interactions in the dataset
for i = 1:length(allMMF)

    % Determine the frame corresponding to copulation time
    if ~isempty(allMMF(i).copTime)
        copFrame = find(allMMF(i).timestamps(1, :) > allMMF(i).copTime, 1);
    else
        copFrame = length(allMMF(i).timestamps);
    end

    % Loop over each male in the interaction
    for male = 1:2
        % Set ID for the other male
        otherID = 3 - male; % Toggles between 1 and 2
        flyCtr = flyCtr + 1;

        % Extract and filter song events
        ss = allMMF(i).singing(male, 1:copFrame) > 0;
        allSongs = regionprops(ss, 'PixelIdxList');
        lengths = arrayfun(@(x) length(x.PixelIdxList), allSongs);
        allSongs(lengths < minSongFrames) = []; % Remove short events
        ss = zeros(size(ss)); % Initialize filtered song signal
        for j = 1:length(allSongs)
            ss(allSongs(j).PixelIdxList) = 1; % Retain valid events
        end

        % Extract and filter flick events
        fs = allMMF(i).flicking(male, 1:copFrame) > 0;
        allFlicks = regionprops(fs, 'PixelIdxList');
        lengths = arrayfun(@(x) length(x.PixelIdxList), allFlicks);
        allFlicks(lengths < minFlickFrames) = []; % Remove short events
        fs = zeros(size(fs)); % Initialize filtered flick signal
        for j = 1:length(allFlicks)
            fs(allFlicks(j).PixelIdxList) = 1; % Retain valid events
        end

        % Positions of self, rival, and female
        selfPos = [allMMF(i).xPos(male, :); allMMF(i).yPos(male, :)];
        rivalPos = [allMMF(i).xPos(otherID, :); allMMF(i).yPos(otherID, :)];
        FPos = [allMMF(i).xPos(3, :); allMMF(i).yPos(3, :)];
        selfPos = selfPos(:, 1:copFrame);
        rivalPos = rivalPos(:, 1:copFrame);
        FPos = FPos(:, 1:copFrame);

        % Compute angles between entities
        relSelfFemAng = computeRelativeAngle(FPos, selfPos, allMMF(i).theta(3, 1:copFrame));
        relRivalFemAng = computeRelativeAngle(FPos, rivalPos, allMMF(i).theta(3, 1:copFrame));
        relFemSelfAng = computeRelativeAngle(selfPos, FPos, allMMF(i).theta(male, 1:copFrame));
        relRivalSelfAng = computeRelativeAngle(selfPos, rivalPos, allMMF(i).theta(male, 1:copFrame));
        relFemRivalAng = computeRelativeAngle(rivalPos, FPos, allMMF(i).theta(otherID, 1:copFrame));
        relSelfRivalAng = computeRelativeAngle(rivalPos, selfPos, allMMF(i).theta(otherID, 1:copFrame));

        % Compute inter-fly distances
        rivalDist = allMMF(i).dist2ani2(1, 1:copFrame);
        femaleDist = allMMF(i).dist2ani3(male, 1:copFrame);
        rivalDist2Fem = allMMF(i).dist2ani3(otherID, 1:copFrame);

        % Assemble raw and filtered predictors
        rawPredictors = [femaleDist', rivalDist', rivalDist2Fem', ...
                         abs(relFemSelfAng)', abs(relRivalSelfAng)', ...
                         abs(relSelfFemAng)', abs(relRivalFemAng)', ...
                         abs(relSelfRivalAng)', abs(relFemRivalAng)', ...
                         ones(length(relRivalFemAng), 1) * flyCtr];
        filteredPredictors = fillmissing(rawPredictors, 'movmean', 5); % Smooth data
        filteredPredictors = filteredPredictors ./ max(filteredPredictors); % Normalize
        filteredPredictors(:, end) = []; % Exclude fly ID column

        % Append predictors and labels to storage
        store.rawPredictors = cat(1, store.rawPredictors, rawPredictors);
        store.filteredPredictors = cat(1, store.filteredPredictors, filteredPredictors);
        store.song = cat(1, store.song, ss');
        store.flick = cat(1, store.flick, fs');
    end
    disp(['Finished pre-processing data for triad #' num2str(i)])
end

% Clean up workspace and prepare for visualization
clearvars -except store*

allPredictors = store.filteredPredictors; 
allPredictors = fillmissing(allPredictors, 'constant', 0);
allPredictorsRaw = store.rawPredictors; 
allPredictorsRaw = fillmissing(allPredictorsRaw, 'constant', 0);
allKeysS = store.song;
allKeysF = store.flick;

% Perform t-SNE dimensionality reduction
disp('Beginning tSNE pipeline...This will take a while. Have a coffee!');
Y = tsne(allPredictors, 'Perplexity', 200, 'Exaggeration', 5, 'Verbose', 1, 'NumPrint', 5);

% Create histogram bins for visualization
plotRes = 0.5; % Adjust granularity of the plot
binsX = round(min(Y(:, 1)) / 10 - 1) * 10:plotRes:round(max(Y(:, 1)) / 10 + 1) * 10;
binsY = round(min(Y(:, 2)) / 10 - 1) * 10:plotRes:round(max(Y(:, 2)) / 10 + 1) * 10;

% Compute 2D bins for the entire dataset
[~, ~, binX] = histcounts(Y(:, 1), binsX);
[~, ~, binY] = histcounts(Y(:, 2), binsY);

% Combine bin indices into a single linear index for each point
validIdx = binX > 0 & binY > 0; % Exclude out-of-bounds points
linearIdx = sub2ind([length(binsX) - 1, length(binsY) - 1], binX(validIdx), binY(validIdx));

% Initialize probability maps as column vectors
tsneMapVec = zeros((length(binsX) - 1) * (length(binsY) - 1), 1);
tsneMapSongVec = zeros((length(binsX) - 1) * (length(binsY) - 1), 1);
tsneMapFlickVec = zeros((length(binsX) - 1) * (length(binsY) - 1), 1);

% Filter valid keys
allKeysS_valid = allKeysS(validIdx);
allKeysF_valid = allKeysF(validIdx);

% Accumulate counts for each map
tsneMapVec = accumarray(linearIdx, double(allKeysS_valid == 0 | allKeysS_valid == 1), size(tsneMapVec));
tsneMapSongVec = accumarray(linearIdx, double(allKeysS_valid == 1), size(tsneMapSongVec));
tsneMapFlickVec = accumarray(linearIdx, double(allKeysF_valid == 1), size(tsneMapFlickVec));

% Reshape vectors back into 2D matrices
tsneMap = reshape(tsneMapVec, [length(binsX) - 1, length(binsY) - 1]);
tsneMapSong = reshape(tsneMapSongVec, [length(binsX) - 1, length(binsY) - 1]);
tsneMapFlick = reshape(tsneMapFlickVec, [length(binsX) - 1, length(binsY) - 1]);

% Normalize and visualize probability maps
pSongMap = fillmissing(tsneMapSong ./ tsneMap, 'constant', 0);
pFlickMap = fillmissing(tsneMapFlick ./ tsneMap, 'constant', 0);

tsneMap = tsneMap ./ sum(tsneMap(:));

% Plot results
subplot(131);
imagesc(imgaussfilt(tsneMap, 1))
title('All positions')
clim([0 .3e-3])

subplot(132);
imagesc(imgaussfilt(pSongMap, 1))
title('P(Song)')
clim([0 .5])

subplot(133);
imagesc(imgaussfilt(pFlickMap, 1))
title('P(Flick)')
clim([0 .5])
colormap gray

%% Get regions of high PSong/Pflick and plot spatial configurations.

% Apply Gaussian filtering to smooth the song probability map. Can change
% to flick. 
allP = imgaussfilt(pFlickMap, 1); % Smoothed song probability map

% Threshold to identify regions of high probability (e.g > 0.3).
bw = allP > 0.4;

% Identify connected regions in the binary mask
conn = regionprops(bw, 'PixelIdxList', 'Area');

% Sort regions by area (descending order) and select the largest 4 regions
areas = arrayfun(@(x) x.Area, conn); % Extract areas
[~, idx] = sort(areas, 'descend');  % Sort areas
conn = conn(idx(1:4));              % Retain the top 4 largest regions

% Initialize binary mask for selected regions
bw = false(size(allP));  % Preallocate logical array
bw(vertcat(conn.PixelIdxList)) = true;  % Logical indexing

% Map raw data points to spatial bins within each region
for i = 1:length(conn)
    % Convert pixel indices to subscripts for rows and columns
    [rows, cols] = ind2sub(size(allP), conn(i).PixelIdxList);
    YVals = binsY(cols); % Y-coordinates of the bins
    XVals = binsX(rows); % X-coordinates of the bins

    % Initialize a container for all raw positions in the current region
    conn(i).allRawPositions = [];
    for j = 1:length(XVals)
        % Select points falling within the current X and Y bins
        XSelect = Y(:,1) > XVals(j) & Y(:,1) < (XVals(j) + plotRes);
        YSelect = Y(:,2) > YVals(j) & Y(:,2) < (YVals(j) + plotRes);
        inBin = find(XSelect & YSelect); % Find indices of points in the bin
        conn(i).allRawPositions = cat(1, conn(i).allRawPositions, inBin); % Aggregate
    end
end

% Reconstruct positions of the three animals (singing male, rival, female)
plotCtr = 1; % Counter for subplot positioning
spatialMaps(length(conn)).m1 = [];  % Preallocate for spatial maps
spatialMaps(length(conn)).rival = [];
for i = 1:length(conn)

    % Extract distance and angle predictors for all animals
    m1_dist2Fem = allPredictorsRaw(conn(i).allRawPositions, 1); % Male1 distance to Female
    rival_dist2Fem = allPredictorsRaw(conn(i).allRawPositions, 3); % Rival distance to Female
    m1_angleFromFem = allPredictorsRaw(conn(i).allRawPositions, 6); % Male1 angle from Female
    rival_angleFromFem = allPredictorsRaw(conn(i).allRawPositions, 7); % Rival angle from Female

    % Calculate centroids of males in the female's coordinate frame
    m1Centroid = [cos(m1_angleFromFem) sin(m1_angleFromFem)] .* m1_dist2Fem;
    rivalCentroid = [cos(rival_angleFromFem) sin(rival_angleFromFem)] .* rival_dist2Fem;
    femaleCentroid = zeros(size(m1Centroid)); % Female fixed at origin

    % Define spatial bin edges
    binsX2 = -200:2.5:200;
    binsY2 = -1:2.5:200;

    % Bin the spatial positions of each male into the defined grid
    [spatialMaps(i).m1, ~, ~] = histcounts2(m1Centroid(:,1), m1Centroid(:,2), binsX2, binsY2);
    [spatialMaps(i).rival, ~, ~] = histcounts2(rivalCentroid(:,1), rivalCentroid(:,2), binsX2, binsY2);

    % Normalize spatial maps to convert counts into probabilities
    spatialMaps(i).m1 = spatialMaps(i).m1 ./ sum(spatialMaps(i).m1(:));
    spatialMaps(i).rival = spatialMaps(i).rival ./ sum(spatialMaps(i).rival(:));

    % Generate colormaps for visualization
    m1Col = [8 106 119] ./ 255; % Color for singing male
    rivalCol = [0 210 179] ./ 255; % Color for rival male
    m1Map = [linspace(1, m1Col(1), 256); linspace(1, m1Col(2), 256); linspace(1, m1Col(3), 256)];
    rivalMap = [linspace(1, rivalCol(1), 256); linspace(1, rivalCol(2), 256); linspace(1, rivalCol(3), 256)];

    % Plot spatial maps
    subplot(length(conn), 2, plotCtr);
    imagesc([fliplr(spatialMaps(i).m1) spatialMaps(i).m1]); % Mirrored heatmap for singing male
    title('Singing Male');
    clim([0 0.0025]);
    colormap(m1Map');
    plotCtr = plotCtr + 1;
    axis off;

    subplot(length(conn), 2, plotCtr);
    imagesc([fliplr(spatialMaps(i).rival) spatialMaps(i).rival]); % Mirrored heatmap for rival male
    title('Rival Male');
    clim([0 0.0025]);
    colormap(rivalMap');
    plotCtr = plotCtr + 1;
    axis off;
end