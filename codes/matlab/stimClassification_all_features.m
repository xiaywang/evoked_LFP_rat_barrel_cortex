%% merged data
startingFolder = pwd;
foldername = uigetdir(startingFolder, 'Select the folder containing all the different data folders');
% layers = ["0240", "0360", "0480", "0720", "0960", "1440", "1680"];
layers = ["0200", "0350", "0500", "0750", "1100", "1500", "1750"];

classnames = [0.5; 1.0; 1.4];

mean_valAcc = []; mean_testAcc = [];
max_valAcc = {}; max_testAcc = {};
max_valAcc_idx = {}; max_testAcc_idx = {};

for l = 1:7
    clear trainData testData time
    [trainData, testData, time] = Preprocessing.merge(foldername, layers(l), 0.8, 0.2);
    % classnames = [0.5; 1; 1.4];

    nTrain = size(trainData, 2);
    nTest = size(testData, 2);

    % preprocessing
    % train data
    sigma = 0.6;
    [trainData, time] = Preprocessing.cropData(trainData, time, 7, 45);
    trainData = Preprocessing.gaussian_smooth(trainData, sigma); % HERE

    % test data
    [testData, time] = Preprocessing.cropData(testData, time, 7, 45);
    testData = Preprocessing.gaussian_smooth(testData, sigma); % HERE

    fileID = fopen('all_features.txt','a');
    fprintf(fileID,'########layer %d########\n', l);
    fclose(fileID);
    disp(strcat('layer', num2str(l)))

    [validationAcc, testAcc] = Utils.all_features(trainData, testData, time, classnames, 'all');
    %     alg_type = 'ensemble';
    %     if layer == 1
    %         n_subm = 2; xcorr = 0; issd = 0; isgrad = 1;
    %     elseif layer == 2
    %         n_subm = 3; xcorr = 0; issd = 0; isgrad = 0;
    %     elseif layer == 3
    %         n_subm = 3; xcorr = 0; issd = 0; isgrad = 0;
    %     elseif layer == 4
    %         n_subm = 3; xcorr = 0; issd = 1; isgrad = 0;
    %     elseif layer == 5
    %         n_subm = 2; xcorr = 0; issd = 1; isgrad = 1;
    %     elseif layer == 6
    %         n_subm = 3; xcorr = 1; issd = 1; isgrad = 1;
    %     elseif layer == 7
    %         n_subm = 2; xcorr = 0; issd = 0; isgrad = 0;
    %     end

    % only with new features???
    % n_subm = 2; xcorr = 0; issd = 0; isgrad = 0;
    %     [validationAcc, testAcc] = Utils.new_features(trainData, testData, classnames, n_subm, xcorr, issd, isgrad, alg_type);

    [max_valAcc{end+1}, max_valAcc_idx{end+1}] = max(validationAcc);
    [max_testAcc{end+1}, max_testAcc_idx{end+1}] = max(testAcc);

    mean_valAcc = cat(2, mean_valAcc, mean(validationAcc));
    mean_testAcc = cat(2, mean_testAcc, mean(testAcc));
end