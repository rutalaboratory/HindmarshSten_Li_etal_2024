%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Basic computation of acoustic features of natural pulse song and
% agonistic song
%
% Related to Cell manuscript "Male-Male Interactions Shape Mate Selection
% in Drosophila" by Tom Hindmarsh Sten, Rufei Li, Florian Hollunder, Shade
% Eleazer, and Vanessa Ruta
%
% Description: 
% This script was provided to demonstrate the parameters and
% selection thresholds we used when computing acoustic featurs of natural
% pulse song and agonistic song. 
% As described in the manuscript, natural courtship pulse song and agonistic 
% song bouts were detected by manually selecting unilateral wing extensions 
% and bilateral wing flicks in MMF assays with one winged and one wingless 
% males. 
% Individual pulses within the song bout were manually detected.
%
% Input: 
% Example structure with pre-processed acoustic recordings (amplitude
% normalized)and their sample rate, courtship pulse song pulses, agonistic song
% pulses
% 
% Output: 
% Inter pulse interval, carrier frequency, pulse number in each
% song bout, relative amplitude comparison between agonistic and courtship
% pulse song
%
%
% Dependencies: 
% Function to get carrier frequency (getPulseFreq.m) from 
% Clemens, Coen, Roemscheid, et al. Current Biology. 2018
% Customized function to calculate distribution: dist.m
%
% Last updated on 2024-12-09 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% compute IPI, pulse train bout length, carrier frequency, and amplitude
%load('Example_natural_song_flicks.mat');

Allfiles = Allfiles_compete;
clearvars -except Allfiles

Ctr = 0;

% set max and min values for inter-pulse-intervals
% set IPI bin values for histogram distribution
IPI_min = 0.01;
IPI_bin = 0.005;
IPI_max = 0.2;

% set max and min values for carrier frequencies
% set carrier frequency bin values for histogram distribution
Fpeak_min = 100;
Fpeak_bin = 10;
Fpeak_max = 1000;

for i = 1:length(Allfiles)

    % only include assays with more than 50 flick bouts for average
    % acoustic statistics comparison
    if length(Allfiles(i).flick) > 50;
        Ctr = Ctr + 1;
    else
        continue;
    end

    % Fs: acoustic recording sample rate
    Fs = Allfiles(i).Fs;

    % acoustic recording of each assay
    y = Allfiles(i).sound;
    % timestamps cooresponding to each frame in acoustic recording
    t = 0:1/Fs:length(Allfiles(i).sound)/Fs;

    % convert pulse event unit from acoustic frame to second
    flickTimes = Allfiles(i).flick/Fs;
    songTimes = Allfiles(i).sing/Fs;

    ind_flick_Ctr = 1;
    ind_sing_Ctr = 1;
    
    for k = 1:2 %k = 1 for flicking & k = 2 for singing
        if k == 1;
            pulse_time = flickTimes; 
            beh_bout = Allfiles(i).Flicks_beh_bouts;
        else
            pulse_time = songTimes;
            beh_bout = Allfiles(i).UWE_beh_bouts;
        end

        threshold = 0.5; % unit: second

        % when pulses are more than 0.5 second apart, they won't be considered as
        % in the same pulse train
    
        % compute IPI based on manually labeled pulse event
        % note: could also use manually labeled behavioral epochs to define
        % pulse train, sometimes, unilateral wing extension epochs may
        % include multiple pulse trains, so still used pulse event to
        % separate trains (lead to similiar results, code commented out below); 
        % also, in this way, the analysis could be consistent with R85A12 activation
        % agonistic song analysis, where we did not label behavioral epochs
        delta = diff(pulse_time);
        indx_pulse = zeros(1,length(pulse_time));
        % select pulses that are not at the start or end of a pulse train
        % pulse trains with less than 3 pulses are naturally excluded in this way
        for j = 2:length(pulse_time)-1
            if delta(j-1) < threshold & delta(j) < threshold
                indx_pulse(j) = 1;
            end
        end

        % IPI
        indx_delta = zeros(1,length(delta));
        indx_delta(find(indx_pulse==1)) = 1;
        indx_delta(find(indx_pulse==1)-1) = 1;
        IPI_ind = delta(find(indx_delta==1));

        % identify start and end of pulse trains to calculate pulse
        % train lengths
        % pulse indices to separate pulse trains
        indices = find(delta > threshold);
        train_start = [1;indices+1];
        train_end = [indices; length(pulse_time)];
        bout_length = train_end - train_start + 1;
        bout_length_mean = mean(bout_length);
        

        %{
        IPI_ind = [];
        bout_length = [];
        for b = 1:length(beh_bout(:,1))
            bout_start = beh_bout(b,1);
            bout_end = beh_bout(b,2);
            pulse_time_indices = find(pulse_time >= bout_start & ...
                pulse_time <= bout_end);

            bout_length(b)=length(pulse_time_indices);

            pulse_time_in_train = pulse_time(pulse_time_indices);
            IPI_bout = diff(pulse_time_in_train);
            IPI_ind = [IPI_ind; IPI_bout];
        end
        bout_length_mean = mean(bout_length);
        %}

        % refine pulse region for carrier frequency and amplitude calculation
        ind_pulse_Ctr = 0;
        ind_pulse = [];
        amplitude_norm = [];
        for n = 1:length(pulse_time)
            %traditional pulse half width to be considered: 0.015 sec
            timeRange = [pulse_time(n)-0.015 pulse_time(n)+0.015];
            idxRange = find(t>timeRange(1)&t<timeRange(end));
            subsound = y(idxRange);

            %compute envelope to extract peak and troughs for pulses
            [ub,lb] = envelope(subsound,5,'peak'); 

            %find the envelope peak with max value in this region, 
            % which indicate the pulse location, and then trace troughs
            [~,peak_all] = findpeaks(ub,'MinPeakProminence',0.03);
            if ~isempty(peak_all)
                [~, sorted_indx] = sort(ub(peak_all));
                center = peak_all(sorted_indx(end));
            end
            % locate troughs
            lms = find(islocalmin(ub));

            if ~isempty(center) && any(lms<center) && any(lms>center) 
                troughs(1) = lms(find(lms<center,1));
                troughs(2) = lms(find(lms>center,1));
               
                actual_pulse = subsound(troughs(1):troughs(2));
                ind_pulse_Ctr = ind_pulse_Ctr + 1;

                % fill the empty array with 0 to align pulses in a matrix
                % note additional 0 in an array won't affect carrier
                % frquency analysis
                empty_pre = 300-(center-troughs(1));
                empty_post = 300-(troughs(2)-center);
                ind_pulse(:,ind_pulse_Ctr) = ...
                    [zeros(empty_pre,1);actual_pulse;zeros(empty_post,1)];
                
                % amplitude
                amplitude_norm(ind_pulse_Ctr) = ub(center);
                % note: each acoustic recording is normalized based on the
                % max amplitude in that recording in the pre-processing
                % step, so we only compare amplitudes within an recording
                % instead of across recordings
               
            end
        end

        % get carrier frequency of all sing/flick pulses in an assay
        % function from Clemens, Coen, Roemscheid, et al. Current Biology.
        % 2018
        [Fpeak_ind,amp,F]=getPulseFreq(ind_pulse,Fs);

        % save output
        if k == 1
            % IPI distribution of flick pulses in each recording
            flick_IPI_prob(Ctr,:) = dist(IPI_ind,IPI_bin,IPI_min,IPI_max);
            % median IPI of flick pulses in each recording
            flick_IPI_median(Ctr) = median(IPI_ind);
            % mean flick pulse train bout lengths in each recording
            flick_bout_length(Ctr) = bout_length_mean;
            % carrier frequencies distribution of flick pulses in each recording
            flick_Fpeak_prob(Ctr,:) = dist(Fpeak_ind,Fpeak_bin,Fpeak_min,Fpeak_max);
            % median carrier frequencies of flick pulses in each recording
            flick_Fpeak_median(Ctr) = median(Fpeak_ind);
            % normalized flick amplitude in this recording
            flick_amp = median(amplitude_norm);
          
        else
            % IPI distribution of song pulses in each recording
            sing_IPI_prob(Ctr,:) = dist(IPI_ind,IPI_bin,IPI_min,IPI_max);
            % median IPI of song pulses in each recording
            sing_IPI_median(Ctr) = median(IPI_ind);
            % mean song pulse train bout lengths in each recording
            sing_bout_length(Ctr) = bout_length_mean;
            % carrier frequencies distribution of song pulses in each recording
            sing_Fpeak_prob(Ctr,:) = dist(Fpeak_ind,Fpeak_bin,Fpeak_min,Fpeak_max);
            % median carrier frequencies of song pulses in each recording
            sing_Fpeak_median(Ctr) = median(Fpeak_ind);
            % normalized song amplitude in this recording
            sing_amp = median(amplitude_norm);

        end
    end
    % amplitude comparison between flick and song pulses in eachrecording
    amplitude_comparison(Ctr) = flick_amp/sing_amp;
end

clearvars -except Allfiles amplitude_comparison... 
    sing_IPI_prob sing_IPI_median sing_bout_length sing_Fpeak_median sing_Fpeak_prob...
    flick_IPI_prob flick_IPI_median flick_bout_length flick_Fpeak_median flick_Fpeak_prob...
    IPI_min IPI_bin IPI_max Fpeak_min Fpeak_bin Fpeak_max

%% plot distribution
flick_color = [54,198,84]./255;
sing_color = [83,174,220]./255;

% IPI distribution
figure('Position',[400,400,600,300])
subplot(1,2,1)
title('IPI distribution');
hold on;
xlim([0,0.2]);
ylim([0,0.6]);

x = IPI_min:IPI_bin:IPI_max-IPI_bin;

standard_error = std(sing_IPI_prob)./sqrt(length(sing_IPI_prob(:,1)));
mean_data = mean(sing_IPI_prob);
plot(x,mean_data,'-','Color', sing_color,'Linewidth', 1);
fill([x,fliplr(x)],[mean_data+standard_error,fliplr(mean_data-standard_error)],...
    sing_color,'FaceAlpha', 0.5,'LineStyle', 'none');

standard_error = std(flick_IPI_prob)./sqrt(length(flick_IPI_prob(:,1)));
mean_data = mean(flick_IPI_prob);
plot(x,mean_data,'-','Color', flick_color,'Linewidth', 1);
fill([x,fliplr(x)],[mean_data+standard_error,fliplr(mean_data-standard_error)],...
    flick_color,'FaceAlpha', 0.5,'LineStyle', 'none');


subplot(1,2,2)
title('Carrier frequency distribution');
hold on
xlim([100,500]);
ylim([0,0.2]);

x = Fpeak_min:Fpeak_bin:Fpeak_max-Fpeak_bin;

standard_error = std(sing_Fpeak_prob)./sqrt(length(sing_Fpeak_prob(:,1)));
mean_data = mean(sing_Fpeak_prob);
plot(x,mean_data,'-','Color', sing_color,'Linewidth', 1);
fill([x,fliplr(x)],[mean_data+standard_error,fliplr(mean_data-standard_error)],...
    sing_color,'FaceAlpha', 0.5,'LineStyle', 'none');

standard_error = std(flick_Fpeak_prob)./sqrt(length(flick_Fpeak_prob(:,1)));
mean_data = mean(flick_Fpeak_prob);
plot(x,mean_data,'-','Color', flick_color,'Linewidth', 1);
fill([x,fliplr(x)],[mean_data+standard_error,fliplr(mean_data-standard_error)],...
    flick_color,'FaceAlpha', 0.5,'LineStyle', 'none');


clear x mean_data standar_error 