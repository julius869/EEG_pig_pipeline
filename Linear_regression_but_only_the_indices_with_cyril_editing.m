% Minimal pipeline for pig data
% All the data have been prepared using the BIDS format
% Cyril Pernet and Julius Søgaard - 2024-

% Add EEGLAB to MATLAB path
eeglab_path = 'C:\Users\juliu\Documents\EEG_pig_project\EEGlab_functions\eeglab2024.2';
addpath(eeglab_path);

% Add FieldTrip to MATLAB path
fieldtrip_path = 'C:\Users\juliu\Documents\EEG_pig_project\FieldTrip_package\fieldtrip-20241211';
addpath(fieldtrip_path);
ft_defaults;  % Initialize FieldTrip

% Close any open figures and clear previous EEGLAB data
close all;
clear variables;
rng('default');

% Variables
rawdata_path = 'C:\Users\juliu\Documents\EEG_pig_project\BIDS_EXPORT_raw'; % Path to BIDS dataset

% Initialize EEGLAB
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

% Check and install required plugins if needed
if ~exist('pop_importbids', 'file')
    plugin_askinstall('bids-matlab-tools', [], 1);
end

% Import BIDS data
[STUDY, EEG] = pop_importbids(rawdata_path, ...
    'outputdir', fullfile(rawdata_path, 'derivatives'), ...
    'bidsevent', 'on', 'bidschanloc', 'on');

STUDY = pop_statparams(STUDY, 'default');
[~,~,AvgChanlocs] = std_prepare_neighbors(STUDY, EEG, 'force', 'on');
channel_info = AvgChanlocs.expected_chanlocs;
save(fullfile(rawdata_path, ['derivatives' filesep 'channel_info.mat']), 'channel_info');

% Plot channel locations
figure('Name', 'Channel locations')
set(gcf, 'Color', 'w', 'InvertHardCopy', 'off', 'units', 'normalized', 'outerposition', [0 0 1 1])
topoplot([], AvgChanlocs.expected_chanlocs, 'style', 'blank', 'electrodes', 'labelpoint', ...
    'chaninfo', AvgChanlocs.expected_chanlocs);
title('EEG channels')
drawnow;

% Load template using FieldTrip
cfg = [];
cfg.dataset = 'C:\Users\juliu\Documents\EEG_pig_project\Banana_template\banana_phantom.edf';
template = ft_preprocessing(cfg);
if any(cellfun(@(x) contains(x,'CZ'),template.label))
    index = find(cellfun(@(x) contains(x,'CZ'),template.label));    
    template.label(index) = [];
    data = template.trial{1}; data(index,:) = [];
    template.trial{1} = data;
end
template_data = template.trial{1}; % Extract template data
fs = template.fsample;            % Sampling frequency


% Processing loop: downsample, clean 50Hz, remove bad channels, interpolate, re-reference, ICA, delete bad segments
for s = 1:size(EEG, 2)
    try
        % Save a copy of the original data before cleaning for visual inspection
        EEG_orig = EEG(s);

        % Handle 'trial_type' field if it exists in EEG(s).event
        if isfield(EEG(s).event, 'trial_type')
            [EEG(s).event.type] = EEG(s).event.trial_type;  % Copy values from trial_type to type
            EEG(s).event = rmfield(EEG(s).event, 'trial_type');  % Remove trial_type field
        end

        if isfield(EEG_orig.event, 'trial_type')
            [EEG_orig.event.type] = EEG_orig.event.trial_type;  % Copy values from trial_type to type
            EEG_orig.event = rmfield(EEG_orig.event, 'trial_type');  % Remove trial_type field
        end

        % Find the index of the 'Cz' channel in the EEG data
        cz_index = find(strcmp({EEG(s).chanlocs.labels}, 'CZ')); % Locate the 'Cz' channel in the chanlocs structure
        disp("Index of 'Cz' channel:");
        disp(cz_index); % Debug: Verify if the index of Cz is being found correctly

        % Check if the 'Cz' channel exists
        if ~isempty(cz_index)
            fprintf('Removing Cz channel for subject %d.\n', s);
            
            % Remove the Cz channel from the EEG data
            EEG(s).data(cz_index, :) = []; % Remove the corresponding row in the data matrix (channel data)
            
            % Remove the Cz channel from the chanlocs structure
            EEG(s).chanlocs(cz_index) = []; % Update the channel labels to exclude 'Cz'
            
            % Update the number of channels to reflect the change
            EEG(s).nbchan = EEG(s).nbchan - 1; % Decrease the number of channels by 1
            
            disp('Updated EEG(s).data size after removing Cz:');
            disp(size(EEG(s).data)); % Debug: Verify the size of EEG data after removal
        else
            fprintf('Cz channel not found for subject %d.\n', s);
        end

        % Display the contents of EEG(s).event
        fprintf('Contents of EEG(s).event:\n');
        for i = 1:length(EEG(s).event)
            fprintf('Event %d:\n', i);
            disp(EEG(s).event(i));  % Display the entire event structure
        end

        % Get lengths of EEG data and template data
        len_eeg = size(EEG(s).data, 2); % Length of EEG data
        len_template = size(template_data, 2); % Length of template data

        % Zero-pad the shorter dataset
        diff_len = abs(len_eeg - len_template);
        if len_eeg > len_template
            % Zero-pad template_data
            template_padded = [template_data, zeros(size(template_data, 1), diff_len)];
        elseif len_template > len_eeg
            % Zero-pad EEG data
            EEG(s).data = [EEG(s).data, zeros(size(EEG(s).data, 1), diff_len)];
        else
            % No padding needed
            template_padded = template_data;
        end

        disp("Length of templated:");
        disp(size(template_padded));

        disp("Length of EEG(s).data:");
        disp(size(EEG(s).data));


        % Compute FFT of the EEG dataset
        eeg_fft = fft(EEG(s).data, [], 2); % FFT of EEG data
        template_fft = fft(template_padded, [], 2); % FFT of padded template


        % Define target frequencies for noise removal
        target_frequencies = [10, 20, 30, 40];


        % Determine the longer length between EEG data and template data
        max_len = max(len_eeg, len_template);
        
        % Compute FFT frequency range using the longest length
        freq_range = (fs / max_len) * (0:floor(max_len / 2 - 1));

        % Retain selected frequencies
        channel_tmp = zeros(size(template_fft));
        disp("size of channel_tmp")
        disp(size(channel_tmp))
        template_fft_filtered = zeros(length(target_frequencies),size(template_fft,2)); % Initialize filtered FFT
        disp("size of template_fft_filtered")
        disp(size(template_fft_filtered))
       
        for f = 1:4
            % Reset channel_tmp at the start of each frequency iteration
            channel_tmp = zeros(size(template_fft)); % Ensure it's cleared for every frequency
            
            for ch = 1:size(template_fft, 1) % Loop through channels
                freq = target_frequencies(f);
                [~, target_index] = min(abs(freq_range - freq)); % Closest freq index
                channel_tmp(ch, target_index) = template_fft(ch, target_index);
                
                % Mirror the values for the two-sided FFT
                if target_index > 1 && target_index < len_eeg / 2
                    mirrored_index = len_eeg - target_index + 2;
                    channel_tmp(ch, mirrored_index) = template_fft(ch, mirrored_index);
                end
            end
            
            % Compute the mean across all channels for this frequency
            template_fft_filtered(f, :) = mean(channel_tmp, 1);
        end
        
        % Exclude intercept term from the regression matrix
        X = template_fft_filtered; % Template FFT without intercept term
        disp("dimensions of X")
        disp(size(X))

        % Perform linear regression with intercept
        B = X' \ eeg_fft'; % Solve regression: B will include B0 (intercept)
        
        % Reconstruct EEG FFT using template and intercept
        fitted_fft = (X' * B)'; % Includes B0 and B2; matches dimensions of eeg_fft_filtered
        
        % Compute residuals
        residuals = abs(eeg_fft) - abs(fitted_fft);
        
        % Transform back to time domain
        cleaned_data = real(ifft(residuals, [], 2)); % Residuals transformed back
        
        % Save the cleaned dataset
        EEG(s).data = cleaned_data; 

        % Plot frequency data for verification
        figure;
        
        % Plot the filtered template spectrum (aligned)
        subplot(3, 1, 1);
        plot(freq_range, abs(template_fft_filtered(:, 1:length(freq_range))));
        title('Filtered Template Spectrum (Aligned)');
        xlabel('Frequency (Hz)');
        ylabel('Amplitude');
        xlim([0 50]);
        

        % Plot the residuals spectrum
        subplot(3, 1, 2);
        plot(freq_range, abs(residuals(1, 1:length(freq_range))));
        title('Residuals Spectrum');
        xlabel('Frequency (Hz)');
        ylabel('Amplitude');
        xlim([0 50]);
        
        % Plot the difference spectrum between EEG and residuals 
        subplot(3, 1, 3);
        difference_spectrum = abs(eeg_fft(1, 1:length(freq_range))) - ...
                              abs(residuals(1, 1:length(freq_range)));
        plot(freq_range, difference_spectrum,'LineWidth',3);
        title('Difference Spectrum Between EEG and residuals');
        xlabel('Frequency (Hz)');
        ylabel('Amplitude Difference');
        xlim([0 50]);


                % Plot the residuals spectrum
        figure('Name', 'Residuals Spectrum');
        for row = 1:size(residuals, 1) % Loop through the rows of the residuals matrix
            subplot(size(residuals, 1), 1, row); % Create one subplot per row
            plot(freq_range, abs(residuals(row, 1:length(freq_range))), 'LineWidth', 1.5);
            title(['Residuals Spectrum - Row ', num2str(row)]);
            xlabel('Frequency (Hz)');
            ylabel('Amplitude');
            xlim([0 50]); % Limit to 0-50 Hz
            grid on;
        end

    catch ME
        fprintf('Error processing dataset %d: %s\n', s, ME.message);
    end
end

