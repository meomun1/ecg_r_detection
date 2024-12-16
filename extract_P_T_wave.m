function [P_wave, T_wave] = extract_P_T_wave(ecg_signal, R_peaks, mean_RR, fs, rangeLimit)
    % Define the windows for P and T waves based on mean RR interval
    P_wave = zeros(1, length(ecg_signal));
    T_wave = zeros(1, length(ecg_signal));

    for i = 1:length(R_peaks)
        % P wave window
        P_start = R_peaks(i) - round(0.2 * mean_RR);
        P_end = R_peaks(i) - round(0.04 * mean_RR);
        if P_start < 1
            P_start = 1;
        end
        if P_end > length(ecg_signal)
            P_end = length(ecg_signal);
        end
        P_wave(P_start:P_end) = ecg_signal(P_start:P_end);

        % T wave window
        T_start = R_peaks(i) + round(0.2 * mean_RR);
        T_end = R_peaks(i) + round(0.4 * mean_RR);
        if T_start < 1
            T_start = 1;
        end
        if T_end > length(ecg_signal)
            T_end = length(ecg_signal);
        end
        T_wave(T_start:T_end) = ecg_signal(T_start:T_end);
    end

    % Normalize the P and T waves
    P_wave = (P_wave - min(P_wave)) / (max(P_wave) - min(P_wave));
    T_wave = (T_wave - min(T_wave)) / (max(T_wave) - min(T_wave));

    % Plot the P and T waves
    figure;

    figure(4); set(4, 'Name', 'P detection');
    subplot(2, 1, 1); plot(P_wave(1:rangeLimit));
    title('\bf7. P wave detection'); ylim([0, 1]);

    figure(5); set(5, 'Name', 'T detection');
    subplot(2, 1, 1); plot(T_wave(1:rangeLimit));
    title('\bf8. T wave detection'); ylim([0, 1]);

end