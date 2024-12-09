%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Computing covariance between component behaviors, within and across
% animals
%
% Related to Cell manuscript "Male-Male Interactions Shape Mate Selection
% in Drosophila" by Tom Hindmarsh Sten, Rufei Li, Florian Hollunder, Shade
% Eleazer, and Vanessa Ruta
%
% % Description:
% This script analyzes the relationship between specific behaviors (e.g., 
% wing flicking, speed, distance to male/female) and spatial interactions in 
% fly pairs. The key objectives include:
%   - Calculating and visualizing cross-correlations between behaviors and
%     distances between flies during interactions.
%   - Investigating how the behaviors of one fly influence those of the other
%     in terms of proximity and activity levels.
%   - Generating time-series plots to explore the dynamic interactions
%     between behaviors and spatial variables.
%   - One SOURCE behavior is defined, and the relation of all other
%     component behaviors to this source behavior (performed at t=0) is
%     computed.
%
% Key Features:
%   - Computes cross-correlation between behaviors such as flicking, speed, 
%     and distance to male/female.
%   - Visualizes the results using time-series plots showing cross-correlation 
%     between behaviors for the "self" and "other" flies.
%
% Inputs:
%   - Behavioral and spatial data from fly interactions, including times 
%     of copulation and distance measurements ("exMMF_CantonS_behavior.mat"
%     contains example data from 5 MMF triads. The loaded variable, allMMF,
%     is a 1x5 structure with 23 fields detailing different behavioral
%     parameters. All behavioral states are classified by JAABA models. For
%     each parameter, there is a 3xNFrames double that tracks each of the
%     three animals. The two males are always rows 1 & 2, while the female
%     is always row 3. 
%
% Outputs:
%   - Time-series plots showing cross-correlations of behaviors for both 
%     flies ("self" and "other").
%
% Key Dependencies:
%   - None specified (general MATLAB functionality).
%
% Time resolution and bin width (in frames) can be adjusted depending on 
% the behavior and the desired analysis window.
%
% Author: Tom Hindmarsh Sten
% Date: 12/2/2024
% Affiliation: The Rockefeller University
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Compute cross-covariances

% Load data file
clear;
load('exMMF_CantonS_behavior.mat')
%%
% Define the source behavior and target behaviors to analyze
sourceBehavior = 'flicking';
targetBehaviors = {'flicking','singing','speed','dist2male','dist2female','optimalPosition', ...
                   'chaseFemale','chaseMale','copulAttempt','approach'};

% Initialize a counter for animals
aniCtr = 1;
binWidth = 300; % Define the bin width for cross-correlation (in frames)

% Loop through each behavioral data entry
for i = 1:length(allMMF)
    
    % Determine the frame corresponding to the copulation time, or the last frame if unavailable
    if ~isempty(allMMF(i).copTime)
        copFrame = find(allMMF(i).timestamps(1,:) > allMMF(i).copTime, 1);
    else
        copFrame = length(allMMF(i).timestamps) - 1; % Use the last frame if no copTime
    end

    % Loop through the two flies (1 and 2)
    for fly = 1:2
        
        % Identify the "other" fly (complement of the current fly)
        otherID = 3 - fly;

        % Loop through each target behavior
        for beh = targetBehaviors

            % Handling 'dist2male' behavior
            if strcmpi(beh{1}, 'dist2male')
                % Select and interpolate missing data for the current fly
                if fly == 1
                    IMD = fillmissing(allMMF(i).dist2ani2(fly, 1:copFrame), 'linear');
                else
                    IMD = fillmissing(allMMF(i).dist2ani1(fly, 1:copFrame), 'linear');
                end
                
                % Cross-correlation of distance to male with the source behavior (flicking)
                self(1).(beh{1})(aniCtr, :) = xcov(IMD, allMMF(i).(sourceBehavior)(fly, 1:copFrame), binWidth, 'normalized');
                other(1).(beh{1})(aniCtr, :) = xcov(IMD, allMMF(i).(sourceBehavior)(fly, 1:copFrame), binWidth, 'normalized');

            % Handling 'dist2female' behavior
            elseif strcmpi(beh{1}, 'dist2female')
                % Interpolate missing data for both self and the other fly
                d2f_self = fillmissing(allMMF(i).dist2ani3(fly, 1:copFrame), 'linear');
                d2f_other = fillmissing(allMMF(i).dist2ani3(otherID, 1:copFrame), 'linear');
                
                % Cross-correlation of distance to female with the source behavior
                self(1).(beh{1})(aniCtr, :) = xcov(d2f_self, allMMF(i).(sourceBehavior)(fly, 1:copFrame) > 0, binWidth, 'normalized');
                other(1).(beh{1})(aniCtr, :) = xcov(d2f_other, allMMF(i).(sourceBehavior)(fly, 1:copFrame) > 0, binWidth, 'normalized');

            % Handling 'speed' behavior
            elseif strcmpi(beh{1}, 'speed')
                % Cross-correlation of speed with the source behavior
                self(1).(beh{1})(aniCtr, :) = xcov(allMMF(i).(beh{1})(fly, 1:copFrame), allMMF(i).(sourceBehavior)(fly, 1:copFrame) > 0, binWidth, 'normalized');
                other(1).(beh{1})(aniCtr, :) = xcov(allMMF(i).(beh{1})(otherID, 1:copFrame), allMMF(i).(sourceBehavior)(fly, 1:copFrame) > 0, binWidth, 'normalized');

            % General case for other behaviors
            else
                % Cross-correlation of the behavior with the source behavior
                self(1).(beh{1})(aniCtr, :) = xcov(allMMF(i).(beh{1})(fly, 1:copFrame) > 0, allMMF(i).(sourceBehavior)(fly, 1:copFrame) > 0, binWidth, 'normalized');
                other(1).(beh{1})(aniCtr, :) = xcov(allMMF(i).(beh{1})(otherID, 1:copFrame) > 0, allMMF(i).(sourceBehavior)(fly, 1:copFrame) > 0, binWidth, 'normalized');
            end
        end
        aniCtr = aniCtr + 1; % Increment animal counter
    end
end

% Extract field names for plotting
f = fieldnames(self);

% Define time vector for plotting (-5 to 5 seconds, sampled at 60Hz)
t = linspace(-5, 5, 601);

% Plot cross-correlation results for the self fly
figure;
for i = 1:length(f)
    subplot(5, 2, i);
    plot(t, self.(f{i})', 'Color', [0.0627 0.5020 0.5882 0.2], 'LineWidth', 0.4); hold on;
    plot(t, nanmean(self.(f{i})), 'Color', [0.0627 0.5020 0.5882], 'LineWidth', 2); hold on;
    plot([0 0], [-0.2 0.2], 'k', 'LineWidth', 0.25); % Plot vertical line at 0
    plot([-3 3], [0 0], 'k', 'LineWidth', 0.25); % Plot horizontal line at 0
    title([f{i} ' (self)']);
    xlabel('Time (sec)');
    ylabel('Cov (norm)');
    xlim([-3 3]); ylim([-0.3 0.3]);
end

% Plot cross-correlation results for the other fly
figure;
for i = 1:length(f)
    subplot(5, 2, i);
    plot(t, other.(f{i})', 'Color', [0.3059 0.6863 0.6353 0.2], 'LineWidth', 0.4); hold on;
    plot(t, nanmean(other.(f{i})), 'Color', [0.3059 0.6863 0.6353], 'LineWidth', 2); hold on;
    plot([0 0], [-0.2 0.2], 'k', 'LineWidth', 0.25); % Plot vertical line at 0
    plot([-3 3], [0 0], 'k', 'LineWidth', 0.25); % Plot horizontal line at 0
    title([f{i} ' (other)']);
    xlabel('Time (sec)');
    ylabel('Cov (norm)');
    xlim([-3 3]); ylim([-0.3 0.3]);
end

% Set figure background color to white
set(gcf, 'Color', 'w');