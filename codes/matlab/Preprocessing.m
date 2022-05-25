classdef Preprocessing
    methods (Static)
        
        function varargout = load_data(foldername, layers, varargin)
            % Load the 16x16x69 matrices, merge them into one single
            % variable and give it as output with time variable.
            % If the percentage of how to divide the data into training set,
            % testing set and validation set is given, then the 2 (training
            % and testing sets) or 3 (plus validation set) output variables
            % will be given according to the percentage of division.
            % The percentages must be a number between 0 and 1, if only one
            % number is given, it is the percentage of training set, if two
            % numbers are given, the sum must be 1 and the first is the 
            % percentage of training set and the second is the percentage
            % of testing set. If three numbers are given, the sum must be 1
            % and the third number is the percentage of validation set.
            % The output data is a 5-D array as follows:
            % (layer, number of samples, x dimension, y dimension, time)
            
            % filenames order:
%             {'short-5ms-stim_amp0.5V-0200um-Group Name #1.mat' }
%             {'short-5ms-stim_amp0.5V-0200um-Group Name #10.mat'}
%             {'short-5ms-stim_amp0.5V-0200um-Group Name #11.mat'}
%             ...
%             {'short-5ms-stim_amp0.5V-0200um-Group Name #19.mat'}
%             {'short-5ms-stim_amp0.5V-0200um-Group Name #2.mat' }
%             {'short-5ms-stim_amp0.5V-0200um-Group Name #20.mat'}
%             {'short-5ms-stim_amp0.5V-0200um-Group Name #21.mat'}
%             ...
%             {'short-5ms-stim_amp0.5V-0200um-Group Name #29.mat'}
%             {'short-5ms-stim_amp0.5V-0200um-Group Name #3.mat' }
%             {'short-5ms-stim_amp0.5V-0200um-Group Name #30.mat'}
%             {'short-5ms-stim_amp0.5V-0200um-Group Name #31.mat'}
%             ...
%             {'short-5ms-stim_amp0.5V-1750um-Group Name #78.mat'}
%             {'short-5ms-stim_amp0.5V-1750um-Group Name #79.mat'}
%             {'short-5ms-stim_amp0.5V-1750um-Group Name #8.mat' }
%             {'short-5ms-stim_amp0.5V-1750um-Group Name #80.mat'}
%             {'short-5ms-stim_amp0.5V-1750um-Group Name #9.mat' }
            
            if length(varargin) == 1
                if varargin{1} > 1
                    error('The division percentage must be between 0 and 1')
                else
                    percTrain = varargin{1};
                    percTest = 1-percTrain;
                    percValidation = 0;
                    disp('no percTest inserted --> percTest = 1-percTrain; no percValidation inserted, default value is zero')
                end
            elseif length(varargin) == 2
                percTrain = varargin{1};
                percTest = varargin{2};
                percValidation = 0;
                if percTrain+percTest ~= 1
                    error('Sum of the percentages must be 1')
                end
                disp('no percValidation inserted, default value is zero')
            elseif length(varargin) == 3
                percTrain = varargin{1};
                percTest = varargin{2};
                percValidation = varargin{3};
            end

            directory = dir(fullfile(foldername,'*.mat'));
            filenames = {directory.name}';
            
            trial_data = [];
            for idx = 1:length(layers)
                layer = layers(idx);
                s=1;
                for i = 1:length(filenames)
                    if strfind(filenames{i}, layer)
                        % NB! the trial names are as follows: #1 #10 #11
                        % #12 #13 ... #2 #3 ...etc...
                        load(strcat(foldername, '\', filenames{i}));
                        trial_data(idx,s,:,:,:) = data;
                        s = s+1;
                    end
                end
            end
            
            nSamples = size(trial_data, 2);
            
            if isempty(varargin)
                varargout{1} = trial_data;
                varargout{2} = time;
            elseif percValidation == 0
                nTrain = round(nSamples * percTrain);
                nTest = nSamples - nTrain;
                nValidation = 0;
                varargout{1} = trial_data(:, 1:nTrain, :, :, :);
                varargout{2} = trial_data(:, nTrain+1:end, :, :, :);
                varargout{3} = time;
            else
                nTrain = round(nSamples * percTrain);
                nTest = round(nSamples * percTest);
                varargout{1} = trial_data(:, 1:nTrain, :, :, :);
                varargout{2} = trial_data(:, nTrain+1:nTrain+nTest, :, :, :);
                varargout{3} = trial_data(:, nTrain+nTest+1:end, :, :, :);
                varargout{4} = time;
            end
            
        end
        
        function varargout = merge(foldername, layer, varargin)
            % merge the data from different stimulation amplitudes
            % output is a 5D array: (stimulation, n_samples, x, y, time)
            
            if ~isstring(layer) || strlength(layer)~=4
                error('Layer must be a string of 4 character')
            end
            disp(strcat('layer ', layer))
            
            d = dir(foldername);
            isub = [d(:).isdir]; %# returns logical vector
            nameFolds = {d(isub).name}';
            nameFolds(ismember(nameFolds,{'.','..'})) = [];
            
            if length(varargin) == 1
                if varargin{1} > 1
                    error('The division percentage must be between 0 and 1')
                else
                    percTrain = varargin{1};
                    percTest = 1-percTrain;
                    percValidation = 0;
                    disp('no percTest inserted --> percTest = 1-percTrain; no percValidation inserted, default value is zero')
                end
            elseif length(varargin) == 2
                percTrain = varargin{1};
                percTest = varargin{2};
                percValidation = 0;
                if percTrain+percTest ~= 1
                    error('Sum of the percentages must be 1')
                end
                disp('no percValidation inserted, default value is zero')
            elseif length(varargin) == 3
                percTrain = varargin{1};
                percTest = varargin{2};
                percValidation = varargin{3};
            end
            
            trial_data = [];
            for fold = 1:length(nameFolds)
                foldername_stim = strcat(foldername, '\', nameFolds{fold}, '\');
                directory = dir(fullfile(foldername_stim,'*.mat'));
                filenames = {directory.name}';
                disp(filenames(1))
                s=1;
                for i = 1:length(filenames)
                    if strfind(filenames{i}, layer)
                        % NB! the trial names are as follows: #1 #10 #11
                        % #12 #13 ... #2 #3 ...etc...
                        load(strcat(foldername_stim, '\', filenames{i}));
                        trial_data(fold, s,:,:,:) = data;
                        s = s+1;
                    end
                end
            end
            
            nSamples = size(trial_data, 2);
            
            if isempty(varargin)
                varargout{1} = trial_data;
                varargout{2} = time;
            elseif percValidation == 0
                nTrain = round(nSamples * percTrain);
                varargout{1} = trial_data(:, 1:nTrain, :, :, :);
                varargout{2} = trial_data(:, nTrain+1:end, :, :, :);
                varargout{3} = time;
            else
                nTrain = round(nSamples * percTrain);
                nTest = round(nSamples * percTest);
                varargout{1} = trial_data(:, 1:nTrain, :, :, :);
                varargout{2} = trial_data(:, nTrain+1:nTrain+nTest, :, :, :);
                varargout{3} = trial_data(:, nTrain+nTest+1:end, :, :, :);
                varargout{4} = time;
            end
            
        end
        
        function [out, time] = cropData(inputData, t, first, last)
            % crop data in time
            % input: input data, first time point, last time point
            % output: cropped data
            nDim = length(size(inputData));
            out(:,:,:,:,:) = inputData(:, :, :, :, first:last);
            if length(t) ~= size(out, nDim)
                time = t(first:last);
            else
                time = t;
            end
        end
        
        function out = gaussian_smooth(inputData, sigma)
            % smooth the image of each frame with gaussian kernel of sigma
            % input: input data, sigma of the gaussian kernel
            % output: smoothed data
            out = zeros(size(inputData));
            for t = 1:size(inputData, 5)
                for l = 1:size(inputData, 1)
                    for s = 1:size(inputData, 2)
                        out(l, s, :, :, t) = imgaussfilt(squeeze(inputData(l, s, :, :, t)), sigma);
                    end
                end
            end  
        end

        function filenames = loadDataNames(foldername)
            directory = dir(fullfile(foldername,'*.mat'));
            filenames = {directory.name}';
        end
        function dataPerLayer(foldername, layers)
            % Load the 16x16x69 matrices, concatenate them per layer and
            % save them in .mat files
            % filenames order:
%             {'short-5ms-stim_amp0.5V-0200um-Group Name #1.mat' }
%             {'short-5ms-stim_amp0.5V-0200um-Group Name #10.mat'}
%             {'short-5ms-stim_amp0.5V-0200um-Group Name #11.mat'}
%             ...
%             {'short-5ms-stim_amp0.5V-0200um-Group Name #19.mat'}
%             {'short-5ms-stim_amp0.5V-0200um-Group Name #2.mat' }
%             {'short-5ms-stim_amp0.5V-0200um-Group Name #20.mat'}
%             {'short-5ms-stim_amp0.5V-0200um-Group Name #21.mat'}
%             ...
%             {'short-5ms-stim_amp0.5V-0200um-Group Name #29.mat'}
%             {'short-5ms-stim_amp0.5V-0200um-Group Name #3.mat' }
%             {'short-5ms-stim_amp0.5V-0200um-Group Name #30.mat'}
%             {'short-5ms-stim_amp0.5V-0200um-Group Name #31.mat'}
%             ...
%             {'short-5ms-stim_amp0.5V-1750um-Group Name #78.mat'}
%             {'short-5ms-stim_amp0.5V-1750um-Group Name #79.mat'}
%             {'short-5ms-stim_amp0.5V-1750um-Group Name #8.mat' }
%             {'short-5ms-stim_amp0.5V-1750um-Group Name #80.mat'}
%             {'short-5ms-stim_amp0.5V-1750um-Group Name #9.mat' }
            
            filenames = Preprocessing.loadDataNames(foldername);
            % layers = ["0240", "0360", "0480", "0720", "0960", "1440", "1680"];
            trial_data = [];
            for idx = 1:length(layers)
                layer = layers(idx);
                disp(layer)
                temp = [];
                for i = 1:length(filenames)
                    if strfind(filenames{i}, layer)
                        % NB! the trial names are as follows: #1 #10 #11
                        % #12 #13 ... #2 #3 ...etc...
                        load(strcat(foldername, '\', filenames{i}));
                        temp = cat(1, temp, data);
                    end
                end
                trial_data(idx,:,:,:) = temp(:,:,:);
                size(trial_data)
            end
            
            clear data

            if ~exist(strcat(foldername, '\dataPerLayer'), 'dir')
                mkdir(strcat(foldername, '\dataPerLayer'))
            end
            for i = 1:length(layers)
                layer = strcat('layer', layers(i));
                data(:,:,:) = trial_data(i,:,:,:);
                save(strcat(foldername, '\dataPerLayer\', strcat(layer, '.mat')), 'data', 'time');
            end
        end
        
        function varargout = load_divide(foldername, layers, percTrain, percTest, percValidation)
            % load the data from each layer and divide them into training
            % set, testing set and validation set. Give back also the time
            % variable.
            % input:
            % - foldername: where the original data are (.mat files)
            % - layers: list of strings, each string must contain 4 digits which is the depth of each recording in um
            % - percTrain: percentage of data to be used for training (default = 0.8)
            % - percTest: percentage of data to be used for testing the
            % trained model (default = 0.2)
            % - percValidation: percentage of validation data (default = 0)
            % output:
            % - data_train, data_test, data_validation: each subset of data concatenated along axis 1
            
            if nargin == 4
                if percTrain+percTest ~= 1
                    error('Sum of the percentages must be 1')
                end
                disp('no percValidation inserted, default value is zero')
                percValidation = 0;
            elseif nargin == 3
                disp('no percTest inserted --> percTest = 1-percTrain; no percValidation inserted, default value is zero')
                percValidation = 0;
                percTest = 1-percTrain;
            elseif nargin == 2
                disp('no percentages inserted, default value are assigned')
                percTrain = 0.8;
                percTest = 0.2;
                percValidation = 0;
            elseif nargin < 2
                error('Insert at least two inputs')
            end
            
            if ~exist(strcat(foldername, '\dataPerLayer'), 'dir')
            	Preprocessing.dataPerLayer(foldername, layers);
            end
            
            folder_dataPerLayer = strcat(foldername, '\dataPerLayer');
            filenames = Preprocessing.loadDataNames(folder_dataPerLayer);
            load(strcat(folder_dataPerLayer, '\', filenames{1}));
            nSamples = size(data, 1)/16;
            clear data
            
            if percValidation == 0
                nTrain = round(nSamples * percTrain);
                nTest = nSamples - nTrain;
                nValidation = 0;
            else
                nTrain = round(nSamples * percTrain);
                nTest = round(nSamples * percTest);
                nValidation = nSamples - nTrain - nTest;
            end
            data_train = []; data_test = []; data_validation = [];
            for i = 1:length(layers)
                layer = layers(i);
                load(strcat(folder_dataPerLayer, '\layer', layer, '.mat'));
                data_train = cat(1, data_train, data(1:16*nTrain, :, :));
                data_test = cat(1, data_test, data(16*nTrain+1:16*(nTrain + nTest), :, :));
                if nValidation ~= 0
                    data_validation = cat(1, data_validation, data(16*(nTrain+nTest):end, :, :));
                end
            end
            rmdir(folder_dataPerLayer, 's')
            
            varargout{1} = data_train; varargout{2} = data_test;
            if isempty(data_validation)
                varargout{3} = time;
            else
                varargout{3} = data_validation;
                varargout{4} = time;
            end
            
        end
        
        function [new_data, ref] = align_per_trial(data, t_start, t_end)
            % align data per trial. Each layer of the same trial is shifted
            % by the same amount. The shift amount is determined by finding
            % the min in signal of all the signals in all layers of each
            % trial.
            
            % find the index in time with lowerst mean for each layer in each trial
            nTrial = size(data, 1)/(16*7);
            tot_time = size(data, 3);
            % mean_trial = zeros(1, nTrial);
            i_trial= zeros(1, nTrial);
            tmp_mean_layer = zeros(1, 7);
            tmp_i_layer = zeros(1,7);
            tmp_mean_t = zeros(1, tot_time);
            for trial = 1:nTrial
                for layer = 1:7
                    down = (16*nTrial*(layer-1)) + 16*trial;
                    up = down - 15;
                    for t = 1:tot_time
                        if t < t_end && t > t_start
                            tmp_mean_t(t) = mean(mean(data(up:down, :, t)));
                        end
                    end
                    [tmp_mean_layer(layer), tmp_i_layer(layer)] = min(tmp_mean_t);
                end
                % [mean_trial(trial), tmp_i] = min(tmp_mean_layer)
                tmp_i = find(tmp_mean_layer == min(tmp_mean_layer));
                i_trial(trial) = tmp_i_layer(tmp_i);
            end

            % shift
            new_data = zeros(size(data));
            ref = round(mean(i_trial));
            for trial = 1:nTrial
                for layer = 1:7
                    up = (16*nTrial*(layer-1)) + (16*trial - 15);
                    down = up + 15;
                    shift = i_trial(trial) - ref;
                    if shift == 0
                        new_data(up:down, :, :) = data(up:down, :, :);
                    elseif shift > 0
                        new_data(up:down, :, 1:end-shift) = data(up:down, :, 1+shift:end);
                    else
                        new_data(up:down, :, 1-shift:end) = data(up:down, :, 1:end+shift);
                    end
                end
            end
        end
        
        function [new_data, ref] = align_per_layer(data, t_start, t_end)
            % align data per layer. Each layer of all the trials is aligned.
            
            % find the index in time with lowerst mean for each layer in each trial
            nTrial = size(data, 1)/(16*7);
            tot_time = size(data, 3);
            % mean_trial = zeros(1, nTrial);
            new_data = zeros(size(data));
            ref = [];
            for layer = 1:7
                tmp_mean_layer = zeros(1, 7);
                tmp_i_layer = zeros(1,7);
                for trial = 1:nTrial
                    down = (16*nTrial*(layer-1)) + 16*trial;
                    up = down - 15;
                    tmp_mean_t = zeros(1, tot_time);
                    for t = 1:tot_time
                        if t < t_end && t > t_start
                            tmp_mean_t(t) = mean(mean(data(up:down, :, t)));
                        end
                    end
                    [tmp_mean_layer(trial), tmp_i_layer(trial)] = min(tmp_mean_t);
                end
                
                % shift
                ref = cat(1, ref, round(mean(tmp_i_layer)));
                for trial = 1:nTrial
                    down = (16*nTrial*(layer-1)) + 16*trial;
                    up = down - 15;
                    shift = tmp_i_layer(trial) - ref(end);
                    if shift == 31
                        disp(trial)
                        disp(layer)
                        disp(tmp_i_layer(trial))
                        disp(ref(end))
                    end
                    if shift == 0
                        new_data(up:down, :, :) = data(up:down, :, :);
                    elseif shift > 0
                        new_data(up:down, :, 1:end-shift) = data(up:down, :, 1+shift:end);
                    else
                        new_data(up:down, :, 1-shift:end) = data(up:down, :, 1:end+shift);
                    end
                end
            end
        end
        
        
        function correct_5th_channel()
            startingFolder = pwd;
            folderName = uigetdir(startingFolder, 'Select a folder containing the data');
            directory = dir(fullfile(folderName,'*.mat'));
            files = {directory.name}';

            for i = 1:length(files)
                disp(strcat('loading and processing the data #', num2str(i)))
                load(strcat(folderName, '\', files{i}));
                for t = 1:size(data,3)
                    A = data(:,:,t);
                    A1 = cat(1, A(1:4,:), A(6:end,:));
                    sd = std2(A1);
                    noise = normrnd(0, sd, [1 16]);
                    for j = 1:16
                        A(5,j) = mean([A(4,j), A(6,j)])+noise(j);
                    end
                    data(:,:,t) = A;
                end
                name1_idx = strfind(files{i}, 'Name.mat');
                if isempty(name1_idx)
                    save(strcat(folderName, '\', files{i}), 'data', 'time')
                else
                    name1 = files{i};
                    save(strcat(folderName, '\', name1(1:name1_idx), 'ame #80.mat'), 'data', 'time')
                    delete(strcat(folderName, '\', files{i}));
                end
            end
        end
        
        function name_number_80()
            startingFolder = pwd;
            folderName = uigetdir(startingFolder, 'Select a folder containing the data');
            directory = dir(fullfile(folderName,'*.mat'));
            files = {directory.name}';

            for i = 1:length(files)
                disp(strcat('loading and processing the data #', num2str(i)))
                name1_idx = strfind(files{i}, 'Name.mat');
                if ~isempty(name1_idx)
                    load(strcat(folderName, '\', files{i}));
                    name1 = files{i};
                    save(strcat(folderName, '\', name1(1:name1_idx), 'ame #80.mat'), 'data', 'time')
                    clear data time
                    delete(strcat(folderName, '\', files{i}));
                end
            end
        end      
    end
end