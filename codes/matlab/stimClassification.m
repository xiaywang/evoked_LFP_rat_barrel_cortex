%% merged data
startingFolder = pwd;
foldername = uigetdir(startingFolder, 'Select the folder containing all the different data folders');
% layers = ["0240", "0360", "0480", "0720", "0960", "1440", "1680"];
% layers = ["0500", "0750", "1750"];
layers = ["0200", "0350", "0500", "0750", "1100", "1500", "1750"];
%layers = ["0500"];
%stimulations = [0.6; 0.8; 1.0; 1.2; 1.4];
stimulations = [0.5; 1.0; 1.4];
[trainData, testData, time] = Preprocessing.merge(foldername, layers(7), 0.8, 0.2);

nTrain = size(trainData, 2); % number of samples for each class
nTest = size(testData, 2); % number of samples for each class


%%
sigma = 0.6;
[trainData, time] = Preprocessing.cropData(trainData, time, 8, 45); %50
trainData = Preprocessing.gaussian_smooth(trainData, sigma); % HERE

[testData, time] = Preprocessing.cropData(testData, time, 8, 45); %50
testData = Preprocessing.gaussian_smooth(testData, sigma); % HERE


%%

RPA = featureExtraction.response_peak_amplitude(trainData, time);
pos_rebound = featureExtraction.positive_rebound(trainData, time);
%ROL = featureExtraction.response_onset_latency(trainData, time);
RPL = featureExtraction.response_peak_latency(trainData, time);
tLFP = featureExtraction.time_norm_LFP(trainData, time);

mean_sd = featureExtraction.mean_sd_grad_each_t(trainData, 'n_submatrices', 2, 'issd', 'off', 'isgrad', 'off');
mean_sd = Utils.clear_zeros_mean_sd_grad(mean_sd);
%xcorr = featureExtraction.xcorr2_all(trainData);
%xcorr = Utils.clear_xcorr(xcorr);
features = cat(2, mean_sd, RPA, pos_rebound, RPL, tLFP);

response = [];

for s = 1:length(stimulations)
    for i = 1:nTrain
        response = cat(1, response, stimulations(s));
    end
end

T_train = cat(2, features, response);
rng(7);

[trainedModel, valAcc] = trainClassifier(T_train, stimulations, 'ensemble', 'SubDim', 44, 'nCycles', 22);

%%

% Predictions on test data

RPA = featureExtraction.response_peak_amplitude(testData, time);
pos_rebound = featureExtraction.positive_rebound(testData, time);
%ROL = featureExtraction.response_onset_latency(testData, time);
RPL = featureExtraction.response_peak_latency(testData, time);
tLFP = featureExtraction.time_norm_LFP(testData, time);

mean_sd = featureExtraction.mean_sd_grad_each_t(testData, 'n_submatrices', 2, 'issd', 'off', 'isgrad', 'off');
mean_sd = Utils.clear_zeros_mean_sd_grad(mean_sd);
%xcorr = featureExtraction.xcorr2_all(testData);
%xcorr = Utils.clear_xcorr(xcorr);

T_test = cat(2, mean_sd, RPA, pos_rebound, RPL, tLFP);

correct = [];
for s = 1:length(stimulations)
    for i = 1:nTest
        correct = cat(1, correct, stimulations(s));
    end
end

%%

yfit_trainedModel = trainedModel.predictFcn(T_test);
% Accuracy
Acc_trainedModel = Utils.accuracy(yfit_trainedModel, correct)

%%
figure
Utils.confusion_matrix(validationPredictions, response, stimulations, 'CV')
% xtickangle(0)
% xticklabels({'0.5V', '1.0V', '1.4V'})
% yticklabels({'0.5V', '1.0V', '1.4V'})
figure
Utils.confusion_matrix(yfit_trainedModel, correct, stimulations, 'Test')
set(gca, 'fontsize', 12)
%% span over the hyper-parameters

% % % subspaceDim_values = 20:2:150;
% % % nCycles = 20:2:50;
% % % valAcc_param = zeros(length(nCycles), length(subspaceDim_values));
% % % testAcc_param = zeros(length(nCycles), length(subspaceDim_values));
% % % nc = 1; % index for cycles, number of cycle
% % % for param1 = nCycles
% % %     sd = 1; % subspace dimension
% % % for param2 = subspaceDim_values
% % % 
% % % 
% % %     %
% % %     rng(7);
% % %     [trainedModel, valAcc] = trainClassifier(T_train, stimulations, 'ensemble', 'SubDim', param2, 'nCycles', param1);
% % %     valAcc_param(nc, sd) = valAcc; % cat(1, valAcc_param, valAcc);
% % %     
% % % %     %%
% % % 
% % % 
% % %     %
% % %     yfit_trainedModel = trainedModel.predictFcn(T_test);
% % %     % Accuracy
% % %     Acc_trainedModel = Utils.accuracy(yfit_trainedModel, correct);
% % % 
% % %     testAcc_param(nc, sd) = Acc_trainedModel; % cat(1, testAcc_param, Acc_trainedModel);
% % % 
% % % 
% % %     sd = sd+1;
% % %     
% % % end
% % % nc = nc+1;
% % % end
% % % 
% % % %%
% % % 
% % % [nc, sd] = Utils.find_max_acc(valAcc_param, testAcc_param)
% % % 
% % % subspaceDim_values(13)
% % % nCycles(2)
% % % 
% % % valAcc_param(2, 13)
% % % testAcc_param(2, 13)

%%
% 
% figure
% plot(nCycles, valAcc_subDim, 'b');
% hold on
% plot(nCycles, testAcc_subDim, 'r');
% 
% %%
% 
% [X,Y] = meshgrid(nCycles,subspaceDim_values);
% Z = valAcc_param';
% surf(X,Y,Z)


%%
% %% old dataset
% RPA_old = featureExtraction.response_peak_amplitude(testData_old_L7);
% mean2_sd2_old = featureExtraction.mean2_sd2_each_t(testData_old_L7, 'n_submatrices', 2, 'issd', 'off', 'isgrad', 'off');
% % xcorr = featureExtraction.xcorr2_all(testData);
% T_test_old = cat(2, mean2_sd2_old, RPA_old);
% 
% %%
% yfit_trainedModel_old = trainedModel.predictFcn(T_test_old);
% % Accuracy
% Acc_trainedModel = sum(yfit_trainedModel_old==1.4*ones(19,1))/19
