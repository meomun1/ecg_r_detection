% Clear our variables
clear ecg samplingrate corrected filtered1 peaks1 filtered2 peaks2 fresult

% Load data sample
plotname = 'Sample 2';
load ecgdemodata2;
whos

% Define the range limit for plotting
rangeLimit = 1000; % Adjust this value as needed

% Remove lower frequencies : 0.5Hz - 40Hz
fresult = fft(ecg); % 16999 for Sample 2
value = round(length(fresult)*5/samplingrate); 
fresult(1 : value) = 0;  % 1-85 for Sample 2
fresult(end - value : end) = 0;  % 16914-16999 for Sample 2
corrected = real(ifft(fresult)); % frequency from 5Hz to 495Hz for Sample 2

% Filter - first pass 
winSize = floor(samplingrate * 571 / 1000); % 571ms
if rem(winSize,2)==0                        
    winSize = winSize+1;
end
filtered1=ecgdemowinmax(corrected, winSize); 

% Scales the filtered signal and applies a threshold filter to detect peaks
peaks1=filtered1/(max(filtered1)/7);
for data = 1:1:length(peaks1)
    if peaks1(data) < 4
        peaks1(data) = 0;
    else
        peaks1(data)=1;
    end
end

% Find minimum distance between two peaks
positions=find(peaks1);
distance=positions(2)-positions(1);
for data=1:1:length(positions)-1
    if positions(data+1)-positions(data)<distance 
        distance=positions(data+1)-positions(data);
    end
end

% Optimizes the filter window size based on the minimum distance between peaks
QRdistance=floor(0.04*samplingrate);
if rem(QRdistance,2)==0
    QRdistance=QRdistance+1;
end
winSize=2*distance-QRdistance;

% Filter - second pass
filtered2=ecgdemowinmax(corrected, winSize);

% Detect Atrial Fibrillation
positionsPeak = find(filtered2 > 0); % Find positions of peaks in filtered2
peakDistances = diff(positionsPeak); % Calculate distances between consecutive peaks
meanDistance = mean(peakDistances); % Calculate the mean distance between peaks
distanceStdDev = std(peakDistances); % Calculate the standard deviation of the distances
stdDevThreshold = 0.1 * meanDistance; % 10% of the mean distance
fprintf('Mean distance between peaks = %.2f\n', meanDistance);
fprintf('Standard deviation of distances = %.2f\n', distanceStdDev);
fprintf('Standard deviation threshold = %.2f\n', stdDevThreshold);

if distanceStdDev < stdDevThreshold
    disp('R detection gives a normal result.');
else
    disp('Fluctuation is too large, possible Atrial Fibrillation.');
end

peaks2=filtered2;
for data=1:1:length(peaks2)
    if peaks2(data)<4
        peaks2(data)=0;
    else
        peaks2(data)=1;
    end
end

% Calculates the average heart rate based on the detected peaks.
positions2=find(peaks2);
distanceBetweenFirstAndLastPeaks = positions2(length(positions2))-positions2(1);
averageDistanceBetweenPeaks = distanceBetweenFirstAndLastPeaks/length(positions2);
averageHeartRate = 60 * samplingrate/averageDistanceBetweenPeaks;
disp('Average Heart Rate = ');
disp(averageHeartRate);


% Create figure - stages of processing
% First figure: Original ECG and FFT Filtered ECG
figure(1); set(1, 'Name', strcat(plotname, ' - Processing Stages 1'));

% Original input ECG data
subplot(2, 1, 1); plot((ecg(1:rangeLimit)-min(ecg(1:rangeLimit)))/(max(ecg(1:rangeLimit))-min(ecg(1:rangeLimit))));
title('\bf1. Original ECG'); ylim([-0.2 1.2]);

% ECG with removed low-frequency component
subplot(2, 1, 2); plot((corrected(1:rangeLimit)-min(corrected(1:rangeLimit)))/(max(corrected(1:rangeLimit))-min(corrected(1:rangeLimit))));
title('\bf2. FFT Filtered ECG'); ylim([-0.2 1.2]);

% Second figure: Filtered ECG - 1st pass and Detected peaks
figure(2); set(2, 'Name', strcat(plotname, ' - Processing Stages 2'));

% Filtered ECG (1-st pass) - filter has default window size
subplot(2, 1, 1); stem((filtered1(1:rangeLimit)-min(filtered1(1:rangeLimit)))/(max(filtered1(1:rangeLimit))-min(filtered1(1:rangeLimit))));
title('\bf3. Filtered ECG - 1^{st} Pass'); ylim([0 1.4]);

% Detected peaks in filtered ECG
subplot(2, 1, 2); stem(peaks1(1:rangeLimit));
title('\bf4. Detected Peaks'); ylim([0 1.4]);

% Third figure: Filtered ECG - 2nd pass and Detected Peaks - Finally
figure(3); set(3, 'Name', strcat(plotname, ' - Processing Stages 3'));

% Filtered ECG (2-d pass) - now filter has optimized window size
subplot(2, 1, 1); stem((filtered2(1:rangeLimit)-min(filtered2(1:rangeLimit)))/(max(filtered2(1:rangeLimit))-min(filtered2(1:rangeLimit))));
title('\bf5. Filtered ECG - 2^d Pass'); ylim([0 1.4]);

% Detected peaks - final result
subplot(2, 1, 2); stem(peaks2(1:rangeLimit));
title('\bf6. Detected Peaks - Finally'); ylim([0 1.4]);

% Fourth figure: P and T waves detection 
[P_wave, T_wave] = extract_P_T_wave(corrected, positions2, meanDistance, samplingrate, rangeLimit);


if rangeLimit > 10000
    % Create figure - result
    figure(5); set(5, 'Name', strcat(plotname, ' - Result'));
    % Plotting ECG in green
    plot((ecg(1:rangeLimit)-min(ecg(1:rangeLimit)))/(max(ecg(1:rangeLimit))-min(ecg(1:rangeLimit))), '-g'); title('\bf Comparative ECG R-Peak Detection Plot');

    % Show peaks in the same picture
    hold on
    % Stemming peaks in dashed black
    stem(peaks2(1:rangeLimit)'.*((ecg(1:rangeLimit)-min(ecg(1:rangeLimit)))/(max(ecg(1:rangeLimit))-min(ecg(1:rangeLimit)))'), ':k');
    % Hold off the figure
    hold off
end
% Detect R-peak
% Bandpass filter (0.5Hz - 50Hz)
filt_low = 0.5; 
filt_high = 40;
[b, a] = butter(2, [filt_low filt_high] / (samplingrate / 2), 'bandpass');
filtered_ecg = filtfilt(b, a, ecg);

% Derivative 
diff_ecg = diff(filtered_ecg);

% Squaring the signal
squared_ecg = diff_ecg .^ 2;

% Moving Window Integration
winSize = floor(0.15 * samplingrate); 
integrated_ecg = conv(squared_ecg, ones(1, winSize) / winSize, 'same');

% Thresholding to detect R peaks
threshold = median(integrated_ecg) + 0.4 * std(integrated_ecg); % Adaptive threshold based on median and std
r_peaks = find(integrated_ecg > threshold);

% Ensure r_peaks is a column vector
r_peaks = r_peaks(:);

% Remove closely spaced peaks (Refractory period handling)
min_dist = 0.2 * samplingrate; 
valid_peaks = r_peaks([true; diff(r_peaks) > min_dist]);

% Refine detected R-peaks to align with actual signal peaks
refined_peaks = [];
search_window = round(0.25 * samplingrate); 
for i = 1:length(valid_peaks)
    % Define a local search window around each detected peak
    peak_range = max(1, valid_peaks(i) - search_window):min(length(filtered_ecg), valid_peaks(i) + search_window);
    
    % Find the maximum point within the filtered ECG signal
    [~, local_max_idx_filtered] = max(filtered_ecg(peak_range));
    
    % Align to the original ECG signal
    peak_range_original = max(1, peak_range(local_max_idx_filtered) - search_window):min(length(ecg), peak_range(local_max_idx_filtered) + search_window);
    [~, local_max_idx_original] = max(ecg(peak_range_original));
    
    % Save the global index of the refined peak
    refined_peaks = [refined_peaks; peak_range_original(local_max_idx_original)];
end
valid_peaks = refined_peaks;

% Display detected R-peaks
disp('Detected R-Peaks:');
disp(valid_peaks);

% Calculate RR intervals 
RR_intervals = diff(valid_peaks) / samplingrate; % in seconds

% Calculate statistics for RR intervals
mean_RR = mean(RR_intervals); 
std_RR = std(RR_intervals);   
rmssd = sqrt(mean(diff(RR_intervals) .^ 2)); 

% Display metrics
fprintf('Mean RR Interval: %.3f s\n', mean_RR);
fprintf('Standard Deviation of RR Intervals: %.3f s\n', std_RR);
fprintf('RMSSD: %.3f s\n', rmssd);

% Define thresholds for rhythm irregularities
regular_threshold = 0.6 * mean_RR; 
irregular_threshold = 0.9;         

% Detect rhythm regularity
if std_RR < regular_threshold
    disp('Rhythm is regular.');
elseif std_RR < irregular_threshold
    disp('Detected Regular Irregularity (e.g., PVC).');
else
    disp('Detected Irregular Irregularity (e.g., Atrial Fibrillation).');
end

% Poincaré plot visualization
figure;
scatter(RR_intervals(1:end-1), RR_intervals(2:end), 'filled');
xlabel('RR Interval n (s)');
ylabel('RR Interval n+1 (s)');
title('Poincaré Plot of RR Intervals');

% Overlay detected R-peaks on original ECG
figure;
plot(ecg, '-b'); hold on;
stem(valid_peaks, ecg(valid_peaks), 'or', 'LineWidth', 1.5);
title('ECG Signal with Detected R-Peaks');
xlabel('Sample');
ylabel('Amplitude');
legend('ECG Signal', 'Detected R-Peaks');
