%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DLC_JAABA_pipeline Step 3
% Classify behavior using JAABA and convert output to a struture for future
% analysis
% Detailed description of JAABA see https://jaaba.sourceforge.net/
%
% Related to Cell manuscript "Male-Male Interactions Shape Mate Selection
% in Drosophila" by Tom Hindmarsh Sten, Rufei Li, Florian Hollunder, Shade
% Eleazer, and Vanessa Ruta
%
% Description: 
% This script was provided along with other scripts to show how we
% converted DeepLabCut pose tracking output to JAABA-classified behavioral 
% epochs.
%
% Input: 
% trx.mat and perframe folder, organized in a way that fits with JAABA
% projects
%
% Output: 
% scores_xxbehavior.mat files and a structure that contains useful
% information for further customized analysis
%
% Dependencies: 
% Setting up JAABA correctly
% JAABA classifiers (provided)
% Parallel Computing Toolbox was used during JAABADetect in our case to 
% speed up the process, but it is not required
%
%
% Last updated on 2024-12-09 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% JAABA detection
addpath(genpath('/Users/rufeili/Documents/MATLAB/JAABA'));

jabfiles = {'/Users/rufeili/Documents/Test/Example_JAABA_classifiers/singing.jab';...
'/Users/rufeili/Documents/Test/Example_JAABA_classifiers/flicking.jab';...
'/Users/rufeili/Documents/Test/Example_JAABA_classifiers/pursueAny.jab';...
'/Users/rufeili/Documents/Test/Example_JAABA_classifiers/optimalPosition.jab';...
'/Users/rufeili/Documents/Test/Example_JAABA_classifiers/copulAttempt.jab';...
'/Users/rufeili/Documents/Test/Example_JAABA_classifiers/chaseMale.jab';...
'/Users/rufeili/Documents/Test/Example_JAABA_classifiers/chaseFemale.jab';...
'/Users/rufeili/Documents/Test/Example_JAABA_classifiers/approach.jab'};

parent = '/Users/rufeili/Documents/Test/Example_MMF';

cd(parent);
children = dir(parent);
children = children([children.isdir] & ~ismember({children.name}, {'.', '..'}));

expdirs = {};

for i = 1:length(children)
    expdirs{i} = [parent,'/',children(i).name];
end

JAABADetect(expdirs,'jabfiles',jabfiles,'forcecompute','True');

%% convert JAABA scores file to a structure (MMF) for future analysis
clear all
parent = '/Users/rufeili/Documents/Test/Example_MMF';

cd(parent);
children = dir(parent);
children = children([children.isdir] & ~ismember({children.name}, {'.', '..'}));

for i = 1:length(children)
    filedir = children(i).name;
    cd(filedir);
    load('trx.mat');
    n_flies = 3;
    
    Allfiles(i).sex(1,:) = 'M';
    Allfiles(i).sex(2,:) = 'M';
    Allfiles(i).sex(3,:) = 'F';
    
    for j = 1:n_flies
        Allfiles(i).xPos(j,:) = trx(j).x;
    end

    for j = 1:n_flies
        Allfiles(i).yPos(j,:) = trx(j).y;
    end

    % all things in unit mm need to be scaled down 10 times given we used
    % the wrong scale factor in all previous steps

    for j = 1:n_flies
        Allfiles(i).x_mm(j,:) = trx(j).x_mm / 10;
    end
    
    for j = 1:n_flies
        Allfiles(i).y_mm(j,:) = trx(j).y_mm / 10;
    end
    
    for j = 1:n_flies
        Allfiles(i).theta(j,:) = trx(j).theta;
    end

    for j = 1:n_flies
        Allfiles(i).abdlen_mm(j,:) = trx(j).abdlen_mm /10;
    end

    for j = 1:n_flies
        Allfiles(i).wingAngL(j,:) = trx(j).wing_anglel;
    end

    for j = 1:n_flies
        Allfiles(i).wingAngR(j,:) = trx(j).wing_angler;
    end

    for j = 1:n_flies
        Allfiles(i).timestamps(j,:) = trx(j).timestamps;
    end
    
    for j = 1:n_flies
        Allfiles(i).speed(j,:) =...
            sqrt(diff(trx(j).x_mm).^2+diff(trx(j).y_mm).^2)./(trx(j).dt*10);
    end
    
    for j = 1:n_flies
        Allfiles(i).dist2ani1(j,:) =...
            sqrt((trx(j).x_mm-trx(1).x_mm).^2+(trx(j).y_mm-trx(1).y_mm).^2)/10;
    end
    
    for j = 1:n_flies
        Allfiles(i).dist2ani2(j,:) =...
            sqrt((trx(j).x_mm-trx(2).x_mm).^2+(trx(j).y_mm-trx(2).y_mm).^2)/10;
    end
    
    for j = 1:n_flies
        Allfiles(i).dist2ani3(j,:) =...
            sqrt((trx(j).x_mm-trx(3).x_mm).^2+(trx(j).y_mm-trx(3).y_mm).^2)/10;
    end
    
    behavior_list = dir('scores_*.mat');
    
    for k = 1:length(behavior_list)
        clear tmp
        load(behavior_list(k).name);
        for j = 1:n_flies
            tmp(j,:)=cell2mat(allScores.scores(j));
        end
        Allfiles = setfield(Allfiles,{i},behaviorName,tmp);
    end
    Allfiles(i).dir = filedir;
    
    cd ..
    
end

% Allfiles(1).manual_copframe = [];
% usually we also manually scored copulation and add it to the structure
% only frames before copulation were used for futher analysis

clearvars -except Allfiles

