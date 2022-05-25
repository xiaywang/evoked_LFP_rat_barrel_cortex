function [trainedClassifier, validationAccuracy] = trainClassifier(trainingData, classnames, varargin)
% [trainedClassifier, validationAccuracy] = trainClassifier(trainingData)
% returns a trained classifier and its accuracy. This code recreates the
% classification model trained in Classification Learner app. Use the
% generated code to automate training the same model with new data, or to
% learn how to programmatically train models.
%
%  Input:
%      trainingData: a matrix with the same number of columns and data type
%       as imported into the app.
%
%  Output:
%      trainedClassifier: a struct containing the trained classifier. The
%       struct contains various fields with information about the trained
%       classifier.
%
%      trainedClassifier.predictFcn: a function to make predictions on new
%       data.
%
%      validationAccuracy: a double containing the accuracy in percent. In
%       the app, the History list displays this overall accuracy score for
%       each model.
%
% Use the code to train the model with new data. To retrain your
% classifier, call the function from the command line with your original
% data or new data as the input argument trainingData.
%
% For example, to retrain a classifier trained with the original data set
% T, enter:
%   [trainedClassifier, validationAccuracy] = trainClassifier(T)
%
% To make predictions with the returned 'trainedClassifier' on new data T2,
% use
%   yfit = trainedClassifier.predictFcn(T2)
%
% T2 must be a matrix containing only the predictor columns used for
% training. For details, enter:
%   trainedClassifier.HowToPredict

% Auto-generated by MATLAB on 04-Dec-2017 12:12:28

s = size(trainingData, 2);
varnames = strings(1, s); predictorNames = strings(1, s-1);
for i = 1:s
    varnames(i) = strcat('column_', num2str(i));
end
predictorNames = varnames(1:end-1);
varnames = cellstr(varnames);
predictorNames = cellstr(predictorNames);

% Extract predictors and response
% This code processes the data into the right shape for training the
% model.
% Convert input to table
inputTable = array2table(trainingData, 'VariableNames', varnames);

predictors = inputTable(:, predictorNames);
response = trainingData(:,end);
isCategoricalPredictor = [false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false, false];

% Train a classifier
% This code specifies all the classifier options and trains the classifier.
if length(classnames) == 7
    ensembleSubspaceDim = 49;
else
    ensembleSubspaceDim = 124; %26 %34
end
nCycles = 30; %66;
for i = 1:length(varargin)
    if strcmp(varargin{i}, 'nCycles')
        nCycles = varargin{i+1};
    elseif strcmp(varargin{i}, 'SubDim')
        ensembleSubspaceDim = varargin{i+1};
    end
end
disp(ensembleSubspaceDim)
disp(nCycles)
if isempty(varargin)
    template = templateSVM(...
        'KernelFunction', 'linear', ...
        'PolynomialOrder', [], ...
        'KernelScale', 'auto', ...
        'BoxConstraint', 1, ...
        'Standardize', true);
    classificationAlg = fitcecoc(...
        predictors, ...
        response, ...
        'Learners', template, ...
        'Coding', 'onevsone', ...
        'ClassNames', classnames);
end
if strcmp(varargin{1}, 'SVM')
    degree = varargin{2};
    if degree == 0
        template = templateSVM(...
            'KernelFunction', 'gaussian', ...
            'PolynomialOrder', [], ...
            'KernelScale', 16, ...
            'BoxConstraint', 1, ...
            'Standardize', true);
    elseif degree == 1
        template = templateSVM(...
            'KernelFunction', 'linear', ...
            'PolynomialOrder', [], ...
            'KernelScale', 'auto', ...
            'BoxConstraint', 1, ...
            'Standardize', true);
    else
        template = templateSVM(...
            'KernelFunction', 'polynomial', ...
            'PolynomialOrder', degree, ...
            'KernelScale', 'auto', ...
            'BoxConstraint', 1, ...
            'Standardize', true);
    
    end
    classificationAlg = fitcecoc(...
        predictors, ...
        response, ...
        'Learners', template, ...
        'Coding', 'onevsone', ...
        'ClassNames', classnames);
elseif strcmp(varargin{1}, 'ensemble')
    subspaceDimension = max(1, min(ensembleSubspaceDim, width(predictors) - 1));
    classificationAlg = fitcensemble(...
        predictors, ...
        response, ...
        'Method', 'Subspace', ...
        'NumLearningCycles', nCycles, ...
        'Learners', 'discriminant', ...
        'NPredToSample', subspaceDimension, ...
        'ClassNames', classnames); %, 'NPrint', 1
end


% Create the result struct with predict function
predictorExtractionFcn = @(x) array2table(x, 'VariableNames', predictorNames);
algPredictFcn = @(x) predict(classificationAlg, x);
trainedClassifier.predictFcn = @(x) algPredictFcn(predictorExtractionFcn(x));

% Add additional fields to the result struct
trainedClassifier.ClassificationAlg = classificationAlg;
trainedClassifier.About = 'This struct is a trained model exported from Classification Learner R2017b.';
trainedClassifier.HowToPredict = sprintf('To make predictions on a new predictor column matrix, X, use: \n  yfit = c.predictFcn(X) \nreplacing ''c'' with the name of the variable that is this struct, e.g. ''trainedModel''. \n \nX must contain exactly 98 columns because this model was trained using 98 predictors. \nX must contain only predictor columns in exactly the same order and format as your training \ndata. Do not include the response column or any columns you did not import into the app. \n \nFor more information, see <a href="matlab:helpview(fullfile(docroot, ''stats'', ''stats.map''), ''appclassification_exportmodeltoworkspace'')">How to predict using an exported model</a>.');


% Perform cross-validation
partitionedModel = crossval(trainedClassifier.ClassificationAlg, 'KFold', 5);

% Compute validation predictions
[validationPredictions, validationScores] = kfoldPredict(partitionedModel);

% Compute validation accuracy
validationAccuracy = 1 - kfoldLoss(partitionedModel, 'LossFun', 'ClassifError');

assignin('base', 'validationPredictions', validationPredictions);
