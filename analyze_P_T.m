function analyze_P_T(P_wave, T_wave, mean_RR, fs, rangeLimit)
    % Sampling interval (time per sample)
    sampling_interval = 1 / fs; % in seconds

    % Initialize variables
    P_properties = struct('Amplitude', [], 'Width', [], 'Asymmetry', []);
    T_properties = struct('Amplitude', [], 'Width', [], 'Asymmetry', []);

    % Normalize the P and T waves for analysis
    P_wave = (P_wave - min(P_wave)) / (max(P_wave) - min(P_wave));
    T_wave = (T_wave - min(T_wave)) / (max(T_wave) - min(T_wave));

    %% Analyze P wave
    P_wave_indices = find(P_wave > 0); % Extract non-zero values
    if ~isempty(P_wave_indices)
        % Biên độ sóng P
        P_amplitude = max(P_wave(P_wave_indices));
        
        % Độ rộng sóng P
        P_width = (length(P_wave_indices) * sampling_interval) * 1000; % Convert to ms
        
        % Đối xứng sóng P
        P_mid_index = floor(length(P_wave_indices) / 2);
        P_left_sum = sum(P_wave(P_wave_indices(1:P_mid_index)));
        P_right_sum = sum(P_wave(P_wave_indices(P_mid_index+1:end)));
        P_asymmetry = abs(P_left_sum - P_right_sum) / (P_left_sum + P_right_sum);
        
        % Ghi nhận các đặc trưng
        P_properties.Amplitude = P_amplitude;
        P_properties.Width = P_width;
        P_properties.Asymmetry = P_asymmetry;
        
        % In kết quả
        fprintf('P Wave Analysis:\n');
        fprintf('Amplitude: %.2f\n', P_amplitude);
        fprintf('Width: %.2f ms\n', P_width);
        if P_asymmetry > 0.2
            fprintf('Asymmetry detected (%.2f). Possible abnormality.\n', P_asymmetry);
        else
            fprintf('Symmetrical shape (%.2f).\n', P_asymmetry);
        end
    else
        fprintf('No P wave detected for analysis.\n');
    end

    %% Analyze T wave
    T_wave_indices = find(T_wave > 0); % Extract non-zero values
    if ~isempty(T_wave_indices)
        % Biên độ sóng T
        T_amplitude = max(T_wave(T_wave_indices));
        
        % Độ rộng sóng T
        T_width = (length(T_wave_indices) * sampling_interval) * 1000; % Convert to ms
        
        % Đối xứng sóng T
        T_mid_index = floor(length(T_wave_indices) / 2);
        T_left_sum = sum(T_wave(T_wave_indices(1:T_mid_index)));
        T_right_sum = sum(T_wave(T_wave_indices(T_mid_index+1:end)));
        T_asymmetry = abs(T_left_sum - T_right_sum) / (T_left_sum + T_right_sum);
        
        % Ghi nhận các đặc trưng
        T_properties.Amplitude = T_amplitude;
        T_properties.Width = T_width;
        T_properties.Asymmetry = T_asymmetry;
        
        % In kết quả
        fprintf('\nT Wave Analysis:\n');
        fprintf('Amplitude: %.2f\n', T_amplitude);
        fprintf('Width: %.2f ms\n', T_width);
        if T_asymmetry > 0.3
            fprintf('Asymmetry detected (%.2f). Possible abnormality.\n', T_asymmetry);
        else
            fprintf('Symmetrical shape (%.2f).\n', T_asymmetry);
        end
    else
        fprintf('No T wave detected for analysis.\n');
    end
end