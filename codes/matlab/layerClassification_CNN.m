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

% preprocessing
[trainData, time] = Preprocessing.cropData(trainData, time, 8, 32);
[testData, time] = Preprocessing.cropData(testData, time, 8, 32);

%% Concatenate the frames



%%

layers = [
    imageInputLayer([16 16 1])
    
    convolution2dLayer(3,8,'Padding','same')
    reluLayer
    
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,16,'Padding','same')
    reluLayer
    
    maxPooling2dLayer(2,'Stride',2)
    
    convolution2dLayer(3,32,'Padding','same')
    reluLayer
    
    fullyConnectedLayer(10)
    softmaxLayer
    classificationLayer];

