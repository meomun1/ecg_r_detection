%   We are processing two data samples to demonstrate two different situations
for demo = 1:2:3 % try 1 and 3
    %   Clear our variables
    clear ecg samplingrate corrected filtered1 peaks1 filtered2 peaks2 fresult

    %   Load data sample
    switch(demo)
        case 1,
            plotname = 'Sample 1';
            load ecgdemodata1;
            whos

            %  demo: size - 1x1, bytes - 8, class - double.
            %  ecg: size - 1x44604, bytes - 356832, class - double.
            %  samplingrate: size - 1x1, bytes - 8, class - double.
            %  plotname: size - 1x1, bytes - 8, class - double.

        case 3,
            plotname = 'Sample 2';
            load ecgdemodata2;
            whos
            % QRdistance: size - 1x1, bytes - 8, class - double.
            % winSize: size - 1x1, bytes - 8, class - double.
            % averageDistanceBetweenPeaks: size - 1x1, bytes - 8, class - double.
            % averageHeartRate: size - 1x1, bytes - 8, class - double.
            % data - size - 1x1, bytes - 8, class - double.
            % demo - size - 1x1, bytes - 8, class - double.
            % distance - size - 1x1, bytes - 8, class - double.
            % distanceBetweenFirstAndLastPeaks - size - 1x1, bytes - 8, class - double.
            % ecg - size - 1x16999, bytes - 135992, class - double.
            % plotname - size - 1x8, bytes - 16, class - char.
            % positions - size - 1x47, bytes - 376, class - double.
            % positions2 - size - 1x47, bytes - 376, class - double.
            % samplingrate - size - 1x1, bytes - 8, class - double.
    end
    
    %   Remove lower frequencies : 0.5Hz - 40Hz by 
    %   Applies a Fast Fourier Transform (FFT) to the ECG data, removes lower frequencies, and then applies an inverse FFT to get the corrected signal
    fresult = fft(ecg); % 44604 and 16999 for Sample 1 and Sample 2 respectively
    value = round(length(fresult)*5/samplingrate); 
    fresult(1 : value) = 0;  % 1-223 for Sample 1 and 1-85 for Sample 2
    fresult(end - value : end) = 0;  % 44382-44604 for Sample 1 and 16914-16999 for Sample 2
    corrected = real(ifft(fresult)); % frequency from 5Hz to 495Hz for Sample 1 and 5Hz to 495Hz for Sample 2
    
    %   Filter - first pass 
    %   The window size is set to 571 ms, which is the default value
    %   Filters the corrected signal using a window size based on the sampling rate
    winSize = floor(samplingrate * 571 / 1000); % 571ms
    if rem(winSize,2)==0                        % If the window size is even, add 1 to make it odd        
        winSize = winSize+1;
    end
    filtered1=ecgdemowinmax(corrected, winSize);%  a windowed maximum filter to detect peaks in the corrected signal
    
    %   Scales the filtered signal and applies a threshold filter to detect peaks
    %   Scale ecg
    peaks1=filtered1/(max(filtered1)/7);
    %   Filter by threshold filter
    for data = 1:1:length(peaks1)
        if peaks1(data) < 4
            peaks1(data) = 0;
        else
            peaks1(data)=1;
        end
    end

    %  Find minimum distance between two peaks
    positions=find(peaks1);
    distance=positions(2)-positions(1);
    %   Returns minimum distance between two peaks
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
    % Filters the corrected signal again using the optimized window size and applies a threshold filter to detect peaks.
    filtered2=ecgdemowinmax(corrected, winSize);

    positionsPeak = find(filtered2 > 0); % Find positions of peaks in filtered2
    peakValues = filtered2(positionsPeak); % Get the values of the peaks
    peakStdDev = std(peakValues); % Calculate the standard deviation of the peak values
    stdDevThreshold = 0.5; % Adjust this threshold based on your requirements
    if peakStdDev < stdDevThreshold
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


    %   Create figure - stages of processing
    figure(demo); set(demo, 'Name', strcat(plotname, ' - Processing Stages'));


    %   Original input ECG data
    subplot(3, 2, 1); plot((ecg-min(ecg))/(max(ecg)-min(ecg)));
    title('\bf1. Original ECG'); ylim([-0.2 1.2]);


    %   ECG with removed low-frequency component
    subplot(3, 2, 2); plot((corrected-min(corrected))/(max(corrected)-min(corrected)));
    title('\bf2. FFT Filtered ECG'); ylim([-0.2 1.2]);


    %   Filtered ECG (1-st pass) - filter has default window size
    subplot(3, 2, 3); stem((filtered1-min(filtered1))/(max(filtered1)-min(filtered1)));
    title('\bf3. Filtered ECG - 1^{st} Pass'); ylim([0 1.4]);


    %   Detected peaks in filtered ECG
    subplot(3, 2, 4); stem(peaks1);
    title('\bf4. Detected Peaks'); ylim([0 1.4]);


    %   Filtered ECG (2-d pass) - now filter has optimized window size
    subplot(3, 2, 5); stem((filtered2-min(filtered2))/(max(filtered2)-min(filtered2)));
    title('\bf5. Filtered ECG - 2^d Pass'); ylim([0 1.4]);


    %   Detected peaks - final result
    subplot(3, 2, 6); stem(peaks2);
    title('\bf6. Detected Peaks - Finally'); ylim([0 1.4]);


    %   Create figure - result
    figure(demo+1); set(demo+1, 'Name', strcat(plotname, ' - Result'));


    %   Plotting ECG in green
    plot((ecg-min(ecg))/(max(ecg)-min(ecg)), '-g'); title('\bf Comparative ECG R-Peak Detection Plot');


    %   Show peaks in the same picture
    hold on
    %   Stemming peaks in dashed black
    stem(peaks2'.*((ecg-min(ecg))/(max(ecg)-min(ecg)))', ':k');
    %   Hold off the figure
    hold off
end