%% layer classification

% load data
startingFolder = pwd;
foldername = uigetdir(startingFolder, 'Select the folder containing the data');
% layers = ["0240", "0360", "0480", "0720", "0960", "1440", "1680"];
layers = ["0200", "0350", "0500", "0750", "1100", "1500", "1750"];
[trainData, testData, time] = Preprocessing.load_data(foldername, layers, 0.8, 0.2);

nTrain = size(trainData, 2);
nTest = size(testData, 2);

classnames = [1; 2; 3; 4; 5; 6; 7];

%% preprocessing

% training set
[trainData, time] = Preprocessing.cropData(trainData, time, 8, 32);
sigma = 0.6;
trainData = Preprocessing.gaussian_smooth(trainData, sigma);

% testing set
[testData, time] = Preprocessing.cropData(testData, time, 8, 32);
testData = Preprocessing.gaussian_smooth(testData, sigma);

%% training sweeping features and algorithms

[validationAcc, testAcc] = Utils.all_features(trainData, testData, time, classnames)
% rng(7) in the function Utils.all_features

%%
max_valAcc = max(validationAcc)
find(max_valAcc == validationAcc)
max_testAcc = max(testAcc)
find(max_testAcc == testAcc)

% 1.4V
% 19 n_subm 2, xcorr 0, issd 1, isgrad 0 -- validationAcc 0.9821, testAcc 0.9821
% 11 n_subm 3, xcorr 0, issd 1, isgrad 0 -- validationAcc 0.9777, testAcc 0.9911
% 58 n_subm 3, xcorr 0, issd 0, isgrad 1 -- validationAcc 0.9955, testAcc 0.9911
% 61 n_subm 3, xcorr 1, issd 0, isgrad 0 -- validationAcc 0.9955, testAcc 0.9911
% 64 n_subm 3, xcorr 1, issd 1, isgrad 1 -- validationAcc 0.9888, testAcc 0.9911
% n_subm 2, xcorr 1, issd 0, isgrad 0 -- validationAcc 0.9888, testAcc 0.9911

% 1.0V
% 19 n_subm 2, xcorr 0, issd 1, isgrad 0 -- validationAcc 0.9821, testAcc 0.9732
% 11 n_subm 3, xcorr 0, issd 1, isgrad 0 -- validationAcc 0.9911, testAcc 0.9821
% 58 n_subm 3, xcorr 0, issd 0, isgrad 1 -- validationAcc 0.9710, testAcc 0.9821
% 61 n_subm 3, xcorr 1, issd 0, isgrad 0 -- validationAcc 0.9777, testAcc 0.9732
% 64 n_subm 3, xcorr 1, issd 1, isgrad 1 -- validationAcc 0.9732, testAcc 0.9643
% n_subm 2, xcorr 1, issd 0, isgrad 0 -- validationAcc 0.9688, testAcc 0.9643

% 0.5V
% 19 n_subm 2, xcorr 0, issd 1, isgrad 0 -- validationAcc 0.9375, testAcc 0.9196
% 11 n_subm 3, xcorr 0, issd 1, isgrad 0 -- validationAcc 0.9196, testAcc 0.9196
% 58 n_subm 3, xcorr 0, issd 0, isgrad 1 -- validationAcc 0.9152, testAcc 0.8929
% 61 n_subm 3, xcorr 1, issd 0, isgrad 0 -- validationAcc 0.9129, testAcc 0.9286
% 64 n_subm 3, xcorr 1, issd 1, isgrad 1 -- validationAcc 0.8929, testAcc 0.9107
% n_subm 2, xcorr 1, issd 0, isgrad 0 -- validationAcc 0.9107, testAcc 0.9107

% old
% 11 n_subm 3, xcorr 0, issd 1, isgrad 0 -- validationAcc 0.7333, testAcc 0.8571
% 19 n_subm 2, xcorr 0, issd 1, isgrad 0 -- validationAcc 0.7048, testAcc 0.8571
% 58 n_subm 3, xcorr 0, issd 0, isgrad 1 -- validationAcc 0.7714, testAcc 0.9286
% 61 n_subm 3, xcorr 1, issd 0, isgrad 0 -- validationAcc 0.8286, testAcc 1.0000
% 64 n_subm 3, xcorr 1, issd 1, isgrad 1 -- validationAcc 0.8667, testAcc 0.9286
% n_subm 2, xcorr 1, issd 0, isgrad 0 -- validationAcc 0.8095, testAcc 0.8929
