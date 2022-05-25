%% layer classification

% load data
startingFolder = pwd;

foldername = uigetdir(startingFolder, 'Select the folder containing the data');
clear startingFolder
% old old dataset
% layers = ["0240", "0360", "0480", "0720", "0960", "1440", "1680"];
% layers = ["0180", "0360", "0540", "0720", "0900", "1440", "1620"];
layers = ["0200", "0350", "0500", "0750", "1100", "1500", "1750"];
[trainData, testData, time] = Preprocessing.load_data(foldername, layers, 0.8, 0.2);

nTrain = size(trainData, 2);
nTest = size(testData, 2);

%% Preprocessing

% training set
[trainData, time] = Preprocessing.cropData(trainData, time, 8, 32);
sigma = 0.6;
trainData = Preprocessing.gaussian_smooth(trainData, sigma);

% testing set
[testData, time] = Preprocessing.cropData(testData, time, 8, 32);
testData = Preprocessing.gaussian_smooth(testData, sigma);

%%
% feature extraction
mean_sd_train = featureExtraction.mean_sd_grad_each_t(trainData, 'n_submatrices', 3, 'issd', 'on', 'isgrad', 'on');
xcorr_train = featureExtraction.xcorr2_all(trainData);
%xcorr_train = [];
features_train = cat(2, mean_sd_train, xcorr_train);
clear mean_sd_train xcorr_train

% correct response
response = [];
for i = 1:7
    response = cat(1, response, i*ones(nTrain,1));
end
%response = {'1', '2', '3', '4', '5a', '5b', '6'}';
T_train = cat(2, features_train, response);

rng(7);
%% Training
% APPS ==> Classification Learner
rng(7)
[trainedModel, valAcc] = trainClassifier(T_train, 1:7, 'ensemble', 'SubDim', 42, 'nCycles', 22)

%%
figure
Utils.confusion_matrix(validationPredictions, response, 1:7, 'CV')

%%
% feature extraction
mean_sd_test = featureExtraction.mean_sd_grad_each_t(testData, 'n_submatrices', 3, 'issd', 'on', 'isgrad', 'on');
xcorr_test = featureExtraction.xcorr2_all(testData);
%xcorr_test = [];
T_test = cat(2, mean_sd_test, xcorr_test);

% correct answers
correct = [];
for i = 1:7
    correct = cat(1, correct, i*ones(nTest,1));
end

%% fit and predict
yfit_trainedModel = trainedModel.predictFcn(T_test); % fitting
Acc_trainedModel = Utils.accuracy(yfit_trainedModel, correct) % accuracy

%%
figure
Utils.confusion_matrix(yfit_trainedModel, correct, 1:7, 'Test')

%% span over the hyper-parameters

subspaceDim_values = 20:2:150;
nCycles = 20:2:50;
valAcc_param = zeros(length(nCycles), length(subspaceDim_values));
testAcc_param = zeros(length(nCycles), length(subspaceDim_values));
nc = 1; % index for cycles, number of cycle
for param1 = nCycles
    sd = 1; % subspace dimension
for param2 = subspaceDim_values


    %
    rng(7);
    [trainedModel, valAcc] = trainClassifier(T_train, 1:7, 'ensemble', 'SubDim', param2, 'nCycles', param1);
    valAcc_param(nc, sd) = valAcc; % cat(1, valAcc_param, valAcc);
    
%     %%


    %
    yfit_trainedModel = trainedModel.predictFcn(T_test);
    % Accuracy
    Acc_trainedModel = Utils.accuracy(yfit_trainedModel, correct);

    testAcc_param(nc, sd) = Acc_trainedModel; % cat(1, testAcc_param, Acc_trainedModel);


    sd = sd+1;
    
end
nc = nc+1;
end

%%

[nc, sd] = Utils.find_max_acc(valAcc_param, testAcc_param)

%%
subspaceDim_values(12)
nCycles(2)

valAcc_param(15, 34)
testAcc_param(15, 34)
