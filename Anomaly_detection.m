% Load the dataset from ecgdemodata1.mat
dataStruct = load('ecgdemodata1.mat');
ecgData = dataStruct.ecg(:); % Extract ECG signal

% For demonstration purposes, assume labels are synthetic
% Create synthetic labels: 1 for normal, 0 for anomalies
labels = ones(size(ecgData)); 
labels(1:100:end) = 0; % Mark some samples as anomalies for testing

% Reshape data for feature extraction if needed
features = ecgData; % ECG data as features

% Normalize the data
minVal = min(features, [], 'all');
maxVal = max(features, [], 'all');
features = (features - minVal) / (maxVal - minVal);

% Split data into training and testing sets
cv = cvpartition(labels, 'HoldOut', 0.2);
trainData = features(training(cv), :);
testData = features(test(cv), :);
trainLabels = labels(training(cv), :);
testLabels = labels(test(cv), :);

% Separate normal and anomalous data
normalTrainData = trainData(trainLabels == 1, :);
anomalousTrainData = trainData(trainLabels == 0, :);
normalTestData = testData(testLabels == 1, :);
anomalousTestData = testData(testLabels == 0, :);

% Define the Autoencoder
inputSize = size(trainData, 2);
hiddenSize1 = 32;
hiddenSize2 = 16;
hiddenSize3 = 8;

encoderLayers = [
    featureInputLayer(inputSize, "Normalization", "none")
    fullyConnectedLayer(hiddenSize1)
    reluLayer
    fullyConnectedLayer(hiddenSize2)
    reluLayer
    fullyConnectedLayer(hiddenSize3)
    reluLayer
];

decoderLayers = [
    fullyConnectedLayer(hiddenSize2)
    reluLayer
    fullyConnectedLayer(hiddenSize1)
    reluLayer
    fullyConnectedLayer(inputSize)
    sigmoidLayer
];

autoencoderLayers = [
    featureInputLayer(inputSize, "Normalization", "none") % Input layer
    fullyConnectedLayer(hiddenSize1)
    reluLayer
    fullyConnectedLayer(hiddenSize2)
    reluLayer
    fullyConnectedLayer(hiddenSize3)
    reluLayer
    fullyConnectedLayer(hiddenSize2) % Start decoder
    reluLayer
    fullyConnectedLayer(hiddenSize1)
    reluLayer
    fullyConnectedLayer(inputSize)
    regressionLayer % Output layer for reconstruction
];


% Train the autoencoder
options = trainingOptions('adam', ...
    'MaxEpochs', 20, ...
    'MiniBatchSize', 512, ...
    'Shuffle', 'every-epoch', ...
    'ValidationData', {testData, testData}, ...
    'Plots', 'training-progress');

autoencoder = trainNetwork(normalTrainData, normalTrainData, autoencoderLayers, options);

% Reconstruction Loss for Training Data
reconstructions = predict(autoencoder, normalTrainData);
trainLoss = mean(abs(reconstructions - normalTrainData), 2);

% Determine Threshold
threshold = mean(trainLoss) + std(trainLoss);

% Visualization of Training Loss
figure;
histogram(trainLoss, 50, 'Normalization', 'pdf');
hold on;
xline(mean(trainLoss), 'g--', 'LineWidth', 2, 'Label', 'Mean');
xline(threshold, 'b--', 'LineWidth', 2, 'Label', 'Threshold');
xlabel('Reconstruction Loss');
ylabel('Density');
title('Training Loss Distribution');
legend('Train Loss', 'Mean', 'Threshold');
grid on;

% Reconstruction Loss for Anomalous Test Data
anomalousReconstructions = predict(autoencoder, anomalousTestData);
anomalousLoss = mean(abs(anomalousReconstructions - anomalousTestData), 2);

% Visualization of Anomalous Loss
figure;
histogram(anomalousLoss, 50, 'Normalization', 'pdf', 'FaceColor', 'r');
hold on;
xline(mean(anomalousLoss), 'g--', 'LineWidth', 2, 'Label', 'Mean');
xline(threshold, 'b--', 'LineWidth', 2, 'Label', 'Threshold');
xlabel('Reconstruction Loss');
ylabel('Density');
title('Anomalous Loss Distribution');
legend('Anomalous Loss', 'Mean', 'Threshold');
grid on;
