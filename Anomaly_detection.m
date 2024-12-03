% Load the ECG data
load('ecgdemodata1.mat');
load('ecgdemodata2.mat');

% Combine both datasets
data = [val1'; val2'];
labels = [ones(size(val1, 2), 1); zeros(size(val2, 2), 1)]; % Assume val1 is normal (1), val2 is abnormal (0)

% Normalize the data
minVal = min(data, [], 'all');
maxVal = max(data, [], 'all');
features = (data - minVal) / (maxVal - minVal);

% Split data into training and testing sets
cv = cvpartition(labels, 'HoldOut', 0.2);
trainData = features(training(cv), :);
testData = features(test(cv), :);
trainLabels = labels(training(cv), :);
testLabels = labels(test(cv), :);

% Separate normal and anomalous data for training
normalTrainData = trainData(trainLabels == 1, :);

% Define the Autoencoder structure
inputSize = size(trainData, 2);
hiddenSize1 = 32;
hiddenSize2 = 16;
hiddenSize3 = 8;

encoderLayers = [
    featureInputLayer(inputSize, "Normalization", "none")
    fullyConnectedLayer(hiddenSize1, "Activation", "relu")
    fullyConnectedLayer(hiddenSize2, "Activation", "relu")
    fullyConnectedLayer(hiddenSize3, "Activation", "relu")
];

decoderLayers = [
    fullyConnectedLayer(hiddenSize2, "Activation", "relu")
    fullyConnectedLayer(hiddenSize1, "Activation", "relu")
    fullyConnectedLayer(inputSize, "Activation", "sigmoid")
];

autoencoderLayers = [
    encoderLayers
    decoderLayers
];

% Training options
options = trainingOptions('adam', ...
    'MaxEpochs', 20, ...
    'MiniBatchSize', 512, ...
    'Shuffle', 'every-epoch', ...
    'ValidationData', {testData, testData}, ...
    'Plots', 'training-progress');

% Train the Autoencoder
autoencoder = trainNetwork(normalTrainData, normalTrainData, autoencoderLayers, options);

% Calculate Reconstruction Loss for training data
trainReconstructions = predict(autoencoder, normalTrainData);
trainLoss = mean(abs(trainReconstructions - normalTrainData), 2);

% Set thresholds for classification
thresholdNormal = mean(trainLoss) + std(trainLoss);
thresholdArrhythmia = mean(trainLoss) + 2 * std(trainLoss);
thresholdAFib = mean(trainLoss) + 3 * std(trainLoss);

% Test data predictions and loss calculation
testReconstructions = predict(autoencoder, testData);
testLoss = mean(abs(testReconstructions - testData), 2);

% Classification based on thresholds
classifications = zeros(size(testLoss));
for i = 1:length(testLoss)
    if testLoss(i) <= thresholdNormal
        classifications(i) = 1; % Normal
    elseif testLoss(i) <= thresholdArrhythmia
        classifications(i) = 2; % Arrhythmia
    elseif testLoss(i) <= thresholdAFib
        classifications(i) = 3; % Atrial Fibrillation
    else
        classifications(i) = 4; % Severe anomaly (AV Block)
    end
end

%  classification 
for i = 1:length(classifications)
    switch classifications(i)
        case 1
            disp(['Sample ' num2str(i) ': Normal']);
        case 2
            disp(['Sample ' num2str(i) ': Arrhythmia']);
        case 3
            disp(['Sample ' num2str(i) ': Atrial Fibrillation']);
        case 4
            disp(['Sample ' num2str(i) ': Severe Anomaly (AV Block)']);
    end
end

% Visualize the distribution of reconstruction loss
figure;
histogram(testLoss, 50, 'Normalization', 'pdf', 'FaceColor', 'r');
hold on;
xline(mean(testLoss), 'g--', 'LineWidth', 2, 'Label', 'Mean');
xline(thresholdNormal, 'b--', 'LineWidth', 2, 'Label', 'Normal Threshold');
xline(thresholdArrhythmia, 'c--', 'LineWidth', 2, 'Label', 'Arrhythmia Threshold');
xline(thresholdAFib, 'm--', 'LineWidth', 2, 'Label', 'AFib Threshold');
xlabel('Reconstruction Loss');
ylabel('Density');
title('Test Loss Distribution and Thresholds');
legend('Test Loss', 'Mean', 'Normal', 'Arrhythmia', 'AFib');
grid on;
 