%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% DLC_JAABA_pipeline Step 2-2
% Convert features derived from pose tracking to a format that fit with
% JAABA projects
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
%
% Input: 
% Features of individual flies and inter-fly parameters derived from DLC
% pose tracking files, format: .csv files
%
% Output: 
% One trx.mat file and perframe folder with abdlen_mm.mat, which were
% organized in a way that fits with JAABA projects
%
%
% Dependencies: 
% None
%
% Last updated on 2024-12-09 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
parent = '/Users/rufeili/Documents/Test/Example_MMF';

children = dir(parent);
children = children([children.isdir] & ~ismember({children.name}, {'.', '..'}));

cd(parent)
for j = 1:length(children)
    cd(children(j).name);

    files = dir;

    Male1File = [];
    Male2File = [];
    FemaleFile = [];
    for k = 1:length(files)

        if contains(files(k).name,'femaleJAABA.xlsx')
            FemaleFile = files(k).name;
        elseif contains(files(k).name,'male1JAABA.xlsx')
            Male1File = files(k).name;
        elseif contains(files(k).name,'male2JAABA.xlsx')
            Male2File = files(k).name;
        end
    end

    if isempty(Male1File) || isempty(FemaleFile) || isempty(Male2File)
        disp('FILE NOT FOUND')
    end

    male1 = readtable(Male1File);
    male1 = male1(:,2:end);
    male2 = readtable(Male2File);
    male2 = male2(:,2:end);
    female = readtable(FemaleFile);
    female = female(:,2:end);
    names = male1.Properties.VariableNames;

    % organize variables to a specified structure trx.mat 
    % which is required by JAABA
    for i = 1:size(male1,2)
        trx(1).(names{i}) = table2array(male1(:,i))';
    end
    for i = 1:size(male2,2)
        trx(2).(names{i}) = table2array(male2(:,i))';
    end
    for i = 1:size(female,2)
        trx(3).(names{i}) = table2array(female(:,i))';
    end

    for i = 1:3
        trx(i).fps = trx(i).fps(1);
        trx(i).off = trx(i).off(1);
        trx(i).endframe = trx(i).endframe(1);
        trx(i).firstframe = trx(i).firstframe(1);
        trx(i).nframes = trx(i).nframes(1);
        trx(i).dt = trx(i).dt(1:end-1);
        trx(i).nwingsdetected = ~isnan(trx(i).wing_anglel)+~isnan(trx(i).wing_angler);

        trx(i).theta = fillmissing(trx(i).theta,'previous');
        trx(i).theta_mm = fillmissing(trx(i).theta_mm,'previous');
        trx(i).x = fillmissing(trx(i).x,'previous');
        trx(i).x_mm = fillmissing(trx(i).x_mm,'previous');
        trx(i).y = fillmissing(trx(i).y,'previous');
        trx(i).y_mm = fillmissing(trx(i).y_mm,'previous');
        trx(i).a = fillmissing(trx(i).a,'previous');
        trx(i).a_mm = fillmissing(trx(i).a_mm,'previous');
        trx(i).b = fillmissing(trx(i).b,'previous');
        trx(i).b_mm = fillmissing(trx(i).b_mm,'previous');
    end
    save('trx','trx')
    
    mkdir('perframe');
    % add feature abdlen_mm, aligning with the JAABA perframe format
    % given that this feature is required for our JAABA classifiers
    data = cell(1,3);
    for s=1:3
        data{s} = trx(s).abdlen_mm(trx(s).firstframe:trx(s).endframe);
    end
    units.numerator = cell(1,0);
    units.denominator = cell(1,0);
    save(fullfile('perframe/','abdlen_mm'),'data','units')
    
    cd ..
end

