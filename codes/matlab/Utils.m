classdef Utils
    methods (Static)
        function p = plot_mean1(data, varargin)
            % plot the average signal in time
            % input: data and all the inputs of matlab built-in function
            % plot
            % output: plot of average signal in time
            mean1 = featureExtraction.mean_sd_grad_each_t(data, 'n_submatrices', 1);
            t = 1:size(data, 3);
            p = plot(t, mean1, varargin{:});
            hold on
            min_value = min(mean1);
            t_min = t(find(mean1==min_value));
            plot(t_min, min_value, 'ob')
%             max_value = max(mean1);
%             t_max = t(find(mean1==max_value));
%             plot(7*ones(1,7), linspace(max_value, min_value, 7), 'm', 'LineWidth', 1)
            
%             minimum = strcat('\leftarrow min value = ', num2str(min_value));
%             text(t_min, min_value, minimum)
%             maximum = strcat('\leftarrow max value = ', num2str(max_value));
%             text(t_max, max_value, maximum, 'HorizontalAlignment','right')
        end
        function plot_onset(x, y)
            % plot the point when the stimulation is given
            % input:
            % - x = time point of the stimulus onset
            % - y = signal value at onset time
            % output: plot a single red point marking the onset of the
            % stimulation
            plot(x, y, 'or', 'MarkerFaceColor', 'r')
%             strAnnotationType = 'arrow';
            
%             set(gcf, 'Units', 'normalized');
            afXAxisLimits = get(gca, 'XLim');
            afYAxisLimits = get(gca, 'YLim');
%             afEndingPoint = [x, y+0.1*(afYAxisLimits(2)-afYAxisLimits(1))];
%             afStartingPoint = [x, y];
%             afAxesDimensionsAndPositions = get(gca, 'Position');
%             fXAxisPosition = afAxesDimensionsAndPositions(1);
%             fYAxisPosition = afAxesDimensionsAndPositions(2);
%             fXAxisLength = afAxesDimensionsAndPositions(3);
%             fYAxisLength = afAxesDimensionsAndPositions(4);
% %             fXonYAxesRatio = fXAxisLength / fYAxisLength;
%             afFigurePosition = get(gcf, 'Position'); % [left bottom width height]
%             fXonYDimensionRatio = afFigurePosition(3) / afFigurePosition(4);
%             afStartingPoint_FU(1) = ( afStartingPoint(1) - afXAxisLimits(1) ) ...
%             / ( afXAxisLimits(2) - afXAxisLimits(1) ) ...
%             * fXAxisLength + fXAxisPosition;
%             afStartingPoint_FU(2) = ( afStartingPoint(2) - afYAxisLimits(1) ) ...
%             / ( afYAxisLimits(2) - afYAxisLimits(1) ) ...
%             * fYAxisLength + fYAxisPosition;
%             afEndingPoint_FU(1) = ( afEndingPoint(1) - afXAxisLimits(1) ) ...
%             / ( afXAxisLimits(2) - afXAxisLimits(1) ) ...
%             * fXAxisLength + fXAxisPosition;
%             afEndingPoint_FU(2) = ( afEndingPoint(2) - afYAxisLimits(1) ) ...
%             / ( afYAxisLimits(2) - afYAxisLimits(1) ) ...
%             * fYAxisLength + fYAxisPosition;
%           % handleToAnnotation1 = 
% %             annotation( ... %             ... handleToFigure, ... uncomment if necessary
% %             strAnnotationType, ...
% %             [afStartingPoint_FU(1) afEndingPoint_FU(1)], ...
% %             [afStartingPoint_FU(2) afEndingPoint_FU(2)]);
%             annotation('textarrow', [afStartingPoint_FU(1) afEndingPoint_FU(1)], ...
%             [afStartingPoint_FU(2) afEndingPoint_FU(2)],'String','stim. ON', 'HeadStyle', 'none')
            str = {'Stim.', 'ON'};
            text(x-0.02*(afXAxisLimits(2)-afXAxisLimits(1)),y-0.05*(afYAxisLimits(2)-afYAxisLimits(1)),str, 'VerticalAlignment', 'top')
        end
        
        function plot_save_stimulations(data, onset_x, onset_y)
            % plot four average signal in time for each stimulation (0.5V,
            % 1.0V and 1.4V) and save them making a subfolder called
            % Plots/stimulations
            % input: data and onset_x, onset_y
            nTrain = size(data, 1)/7/16;
            for layer = 0:6
                h = figure;
                subplot 311
                top05 = 16*nTrain*layer;
                hold on
                t1 = Utils.plot_mean1(data(top05+1:top05+16, :, :), '.r');
                t2 = Utils.plot_mean1(data(top05+2*16:top05+2*16+16, :, :), '-.k');
                t3 = Utils.plot_mean1(data(top05+3*16:top05+3*16+16, :, :), '--b');
                t4 = Utils.plot_mean1(data(top05+4*16:top05+4*16+16, :, :), ':c', 'LineWidth', 1);
                Utils.plot_onset(onset_x, onset_y)
                legend([t1, t2, t3, t4], {'trial 1', 'trial 2', 'trial 3', 'trial 4'});
                ylabel('Voltage (V)')
                title('0.5V stimulation')
                subplot 312
                top10 = top05+nTrain*16/3;
                hold on
                t1 = Utils.plot_mean1(data(top10+1:top10+16, :, :), '.r');
                t2 = Utils.plot_mean1(data(top10+2*16:top10+2*16+16, :, :), '-.k');
                t3 = Utils.plot_mean1(data(top10+3*16:top10+3*16+16, :, :), '--b');
                t4 = Utils.plot_mean1(data(top10+4*16:top10+4*16+16, :, :), ':c', 'LineWidth', 1);
                Utils.plot_onset(onset_x, onset_y)
                legend([t1, t2, t3, t4], {'trial 1', 'trial 2', 'trial 3', 'trial 4'});
                ylabel('Voltage (V)')
                title('1.0V stimulation')
                subplot 313
                top14 = top10+nTrain*16/3;
                hold on
                t1 = Utils.plot_mean1(data(top14+1:top14+16, :, :), '.r');
                t2 = Utils.plot_mean1(data(top14+2*16:top14+2*16+16, :, :), '-.k');
                t3 = Utils.plot_mean1(data(top14+3*16:top14+3*16+16, :, :), '--b');
                t4 = Utils.plot_mean1(data(top14+4*16:top14+4*16+16, :, :), ':c', 'LineWidth', 1);
                plot(onset_x, onset_y, 'or', 'MarkerFaceColor', 'r')
                Utils.plot_onset(onset_x, onset_y)
                legend([t1, t2, t3, t4], {'trial 1', 'trial 2', 'trial 3', 'trial 4'});
                ylabel('Voltage (V)')
                xlabel('Time points (-)')
                title('1.4V stimulation')
            %     currentFigure = gcf;
            %     title(currentFigure.Children(end), strcat('Layer', num2str(layer+1)));
                ha = axes('Position',[0 0 1 1],'Xlim',[0 1],'Ylim',[0 
                1],'Box','off','Visible','off','Units','normalized', 'clipping' , 'off');

                text(0.5, 1, strcat('Layer',' ', num2str(layer+1)),'HorizontalAlignment' ...
                ,'center','VerticalAlignment', 'top', 'fontsize', 14);

                if ~exist([pwd '\Plots\stimulations'], 'dir')
                    mkdir([pwd '\Plots\stimulations'])
                end
            %     savefig(h, [pwd strcat('/Plots/stimulations/layer', num2str(layer+1), '.fig')])
            %     fig = gcf;
            %     fig.PaperPositionMode = 'auto';
            %     print('-bestfit',[pwd strcat('/Plots/stimulations/bestfit_layer', num2str(layer+1))],'-dpdf')
            %     print('-fillpage','FillPageFigure','-dpdf')
                fig = gcf;
                fig.PaperUnits = 'inches';
                fig.PaperPosition = [0 0 7 8];
                print([pwd strcat('/Plots/stimulations/5by3_layer', num2str(layer+1))],'-dpng','-r0')
            end
        end
        
        function acc = accuracy(yfit, correct)
            % calculate the accuracy of the prediction taking as inputs the
            % fitted values and the correct responses
            result = yfit == correct;
            acc = sum(result)/length(result);
        end
        
        function [validationAcc, testAcc] = all_features(trainData, testData, time, classnames, varargin)
            nTrain = size(trainData, 2);
            nTest = size(testData, 2);
            if classnames(1) == 1
                % correct response
                response = []; % for cross-validation
                for i = 1:7
                    response = cat(1, response, i*ones(nTrain,1));
                end
                % correct answers
                correct = []; % for testing
                for i = 1:7
                    correct = cat(1, correct, i*ones(nTest,1));
                end
            else
                response = [];
                for s = 1:length(classnames)
                    for i = 1:nTrain
                        response = cat(1, response, classnames(s));
                    end
                end
                correct = [];
                for s = 1:length(classnames)
                    for i = 1:nTest
                        correct = cat(1, correct, classnames(s));
                    end
                end
            end
            %time = evalin('base', 'time');
            validationAcc = []; testAcc = [];
            off_on = {'off'; 'on'};
            for alg = 0:4
                fileID = fopen('all_features.txt','a');
                if alg == 0
                    alg_type = 'SVM';
                    fprintf(fileID,'SVM medium gaussian\n');
                elseif alg <= 3
                    alg_type = 'SVM';
                    fprintf(fileID,'SVM degree %d\n', alg);
                else
                    alg_type = 'ensemble'
                    fprintf(fileID,'ensemble\n');
                end
                fclose(fileID);
                for n_subm = 1:3
                    disp(strcat('n_submatrices = ', num2str(n_subm)))
                    for i_xcorr = 0:1
                        disp(strcat('xcorr', num2str(i_xcorr)))
                        for i_sd = 1:2
                            for i_grad = 1:2
                                disp(strcat('i_sd = ', off_on{i_sd}))
                                disp(strcat('i_grad = ', off_on{i_grad}))
                                % feature extraction
                                if classnames(1) == 1 || strcmp(varargin{1}, 'old') %isempty(varargin)
                                    mean_sd_train = featureExtraction.mean_sd_grad_each_t(trainData, 'n_submatrices', n_subm, 'issd', off_on{i_sd}, 'isgrad', off_on{i_grad});
                                    mean_sd_train = Utils.clear_zeros_mean_sd_grad(mean_sd_train);
                                    if i_xcorr == 0
                                        xcorr_train = [];
                                    else
                                        xcorr_train = featureExtraction.xcorr2_all(trainData);
                                        xcorr_train = Utils.clear_xcorr(xcorr_train);
                                    end

                                    features_train = cat(2, mean_sd_train, xcorr_train);

                                    T_train = cat(2, features_train, response);

                                    rng(7)
                                    [trainedModel, valAcc] = trainClassifier(T_train, classnames, alg_type, alg);

                                    validationAcc = cat(2, validationAcc, valAcc);
                                    
                                    mean_sd_test = featureExtraction.mean_sd_grad_each_t(testData, 'n_submatrices', n_subm, 'issd', off_on{i_sd}, 'isgrad', off_on{i_grad});
                                    mean_sd_test = Utils.clear_zeros_mean_sd_grad(mean_sd_test);
                                    if i_xcorr == 0
                                        xcorr_test = [];
                                    else
                                        xcorr_test = featureExtraction.xcorr2_all(testData);
                                        xcorr_test = Utils.clear_xcorr(xcorr_test); 
                                    end

                                    T_test = cat(2, mean_sd_test, xcorr_test);

                                    % fit and predict
                                    yfit_trainedModel = trainedModel.predictFcn(T_test); % fitting
                                    Acc_trainedModel = Utils.accuracy(yfit_trainedModel, correct); % accuracy

                                    testAcc= cat(2, testAcc, Acc_trainedModel);
                                    
                                    fileID = fopen('all_features.txt','a');
                                    fprintf(fileID,'n_subm %d, xcorr %d, issd %d, isgrad %d -- validationAcc %4.4f, testAcc %4.4f\n', n_subm, i_xcorr, i_sd-1, i_grad-1, valAcc, Acc_trainedModel);
                                    fclose(fileID);

                                elseif strcmp(varargin{1}, 'all')
                                    
                                    RPA_train_tmp = featureExtraction.response_peak_amplitude(trainData, time);
                                    RPA_test_tmp = featureExtraction.response_peak_amplitude(testData, time);
                                    pos_reb_train_tmp = featureExtraction.positive_rebound(trainData, time);
                                    pos_reb_test_tmp = featureExtraction.positive_rebound(testData, time);
                                    ROL_train_tmp = featureExtraction.response_onset_latency(trainData, time);
                                    ROL_test_tmp = featureExtraction.response_onset_latency(testData, time);
                                    RPL_train_tmp = featureExtraction.response_peak_latency(trainData, time);
                                    RPL_test_tmp = featureExtraction.response_peak_latency(testData, time);
                                    tLFP_train_tmp = featureExtraction.time_norm_LFP(trainData, time);
                                    tLFP_test_tmp = featureExtraction.time_norm_LFP(testData, time);
                                    
                                    for i_RPA = 0:1
                                        for i_pos_reb = 0:1
                                            for i_ROL = 0:1
                                                for i_RPL = 0:1
                                                    for i_tLFP = 0:1
%                                                         if i_RPA == 0 && i_pos_reb == 0 && i_ROL == 0 && i_RPL == 0 && i_tLFP == 0
%                                                             break
%                                                         end
                                                        disp('new features\n')
                                                        if i_RPA == 0
                                                            RPA_train = [];
                                                            RPA_test = [];
                                                        else
                                                            RPA_train = RPA_train_tmp;
                                                            RPA_test = RPA_test_tmp;
                                                        end
                                                        if i_pos_reb == 0
                                                            pos_reb_train = [];
                                                            pos_reb_test = [];
                                                        else
                                                            pos_reb_train = pos_reb_train_tmp;
                                                            pos_reb_test = pos_reb_test_tmp;
                                                        end
                                                        if i_ROL == 0
                                                            ROL_train = [];
                                                            ROL_test = [];
                                                        else
                                                            ROL_train = ROL_train_tmp;
                                                            ROL_test = ROL_test_tmp;
                                                        end
                                                        if i_RPL == 0
                                                            RPL_train = [];
                                                            RPL_test = [];
                                                        else
                                                            RPL_train = RPL_train_tmp;
                                                            RPL_test = RPL_test_tmp;
                                                        end
                                                        if i_tLFP == 0
                                                            tLFP_train = [];
                                                            tLFP_test = [];
                                                        else
                                                            tLFP_train = tLFP_train_tmp;
                                                            tLFP_test = tLFP_test_tmp;
                                                        end
                                                        
                                                        mean_sd_train = featureExtraction.mean_sd_grad_each_t(trainData, 'n_submatrices', n_subm, 'issd', off_on{i_sd}, 'isgrad', off_on{i_grad});
                                                        mean_sd_train = Utils.clear_zeros_mean_sd_grad(mean_sd_train);
                                                        if i_xcorr == 0
                                                            xcorr_train = [];
                                                        else
                                                            xcorr_train = featureExtraction.xcorr2_all(trainData);
                                                            xcorr_train = Utils.clear_xcorr(xcorr_train); 
                                                        end

                                                        features_train = cat(2, mean_sd_train, xcorr_train, RPA_train, pos_reb_train, ROL_train, RPL_train, tLFP_train);

                                                        T_train = cat(2, features_train, response);

                                                        rng(7)
                                                        [trainedModel, valAcc] = trainClassifier(T_train, classnames, alg_type, alg);

                                                        validationAcc = cat(2, validationAcc, valAcc);

                                                        mean_sd_test = featureExtraction.mean_sd_grad_each_t(testData, 'n_submatrices', n_subm, 'issd', off_on{i_sd}, 'isgrad', off_on{i_grad});
                                                        mean_sd_test = Utils.clear_zeros_mean_sd_grad(mean_sd_test);
                                                        if i_xcorr == 0
                                                            xcorr_test = [];
                                                        else
                                                            xcorr_test = featureExtraction.xcorr2_all(testData);
                                                            xcorr_test = Utils.clear_xcorr(xcorr_test); 
                                                        end

                                                        T_test = cat(2, mean_sd_test, xcorr_test, RPA_test, pos_reb_test, ROL_test, RPL_test, tLFP_test);



                                                        % fit and predict
                                                        yfit_trainedModel = trainedModel.predictFcn(T_test); % fitting
                                                        Acc_trainedModel = Utils.accuracy(yfit_trainedModel, correct); % accuracy

                                                        testAcc= cat(2, testAcc, Acc_trainedModel);
                                                        fileID = fopen('all_features.txt','a');
                                                        fprintf(fileID,'n_subm %d, xcorr %d, issd %d, isgrad %d\n', n_subm, i_xcorr, i_sd-1, i_grad-1);
                                                        fprintf(fileID,'RPA %d, pos_reb %d, ROL %d, RPL %d, tLFP %d -- validationAcc %4.4f, testAcc %4.4f\n', i_RPA, i_pos_reb, i_ROL, i_RPL, i_tLFP,valAcc, Acc_trainedModel);
                                                        fclose(fileID);
                                                    end
                                                end
                                            end
                                        end
                                    end
                                end
                                
                            end
                        end
                    end
                end
            end
        end
        
        function [validationAcc, testAcc] = new_features(trainData, testData, classnames, n_subm, i_xcorr, i_sd, i_grad, alg_type)
            alg = 0;
            if classnames(1) == 1
                nTrain = size(trainData, 1)/7/16;
                nTest = size(testData, 1)/7/16;
                % correct response
                response = [];
                for i = 1:7
                    response = cat(1, response, i*ones(nTrain,1));
                end
                % correct answers
                correct = [];
                for i = 1:7
                    correct = cat(1, correct, i*ones(nTest,1));
                end
            else
                nTrain = size(trainData, 1)/16;
                nTest = size(testData, 1)/16;
                response = [];
                for s = 1:length(classnames)
                    for i = 1:nTrain/length(classnames)
                        response = cat(1, response, classnames(s));
                    end
                end
                correct = [];
                for s = 1:length(classnames)
                    for i = 1:nTest/length(classnames)
                        correct = cat(1, correct, classnames(s));
                    end
                end
            end
            time = evalin('base', 'time');
            validationAcc = []; testAcc = [];
            off_on = {'off'; 'on'};
            
            RPA_train_tmp = featureExtraction.response_peak_amplitude(trainData);
            RPA_test_tmp = featureExtraction.response_peak_amplitude(testData);
            pos_reb_train_tmp = featureExtraction.positive_rebound(trainData);
            pos_reb_test_tmp = featureExtraction.positive_rebound(testData);
            ROL_train_tmp = featureExtraction.response_onset_latency(trainData, time);
            ROL_test_tmp = featureExtraction.response_onset_latency(testData, time);
            RPL_train_tmp = featureExtraction.response_peak_latency(trainData, time);
            RPL_test_tmp = featureExtraction.response_peak_latency(testData, time);
            tLFP_train_tmp = featureExtraction.time_norm_LFP(trainData, time);
            tLFP_test_tmp = featureExtraction.time_norm_LFP(testData, time);
            
            for i_RPA = 0:1
                for i_pos_reb = 0:1
                    for i_ROL = 0:1
                        for i_RPL = 0:1
                            for i_tLFP = 0:1
                                if i_RPA == 0 && i_pos_reb == 0 && i_ROL == 0 && i_RPL == 0 && i_tLFP == 0
                                    continue
                                end
                                disp('new features\n')
                                if i_RPA == 0
                                    RPA_train = [];
                                    RPA_test = [];
                                else
                                    RPA_train = RPA_train_tmp;
                                    RPA_test = RPA_test_tmp;
                                end
                                if i_pos_reb == 0
                                    pos_reb_train = [];
                                    pos_reb_test = [];
                                else
                                    pos_reb_train = pos_reb_train_tmp;
                                    pos_reb_test = pos_reb_test_tmp;
                                end
                                if i_ROL == 0
                                    ROL_train = [];
                                    ROL_test = [];
                                else
                                    ROL_train = ROL_train_tmp;
                                    ROL_test = ROL_test_tmp;
                                end
                                if i_RPL == 0
                                    RPL_train = [];
                                    RPL_test = [];
                                else
                                    RPL_train = RPL_train_tmp;
                                    RPL_test = RPL_test_tmp;
                                end
                                if i_tLFP == 0
                                    tLFP_train = [];
                                    tLFP_test = [];
                                else
                                    tLFP_train = tLFP_train_tmp;
                                    tLFP_test = tLFP_test_tmp;
                                end

                                mean_sd_train = featureExtraction.mean_sd_grad_each_t(trainData, 'n_submatrices', n_subm, 'issd', off_on{i_sd+1}, 'isgrad', off_on{i_grad+1});
                                mean_sd_train = Utils.clear_zeros_mean_sd_grad(mean_sd_train);
                                if i_xcorr == 0
                                    xcorr_train = [];
                                else
                                    xcorr_train = featureExtraction.xcorr2_all(trainData);
                                    xcorr_train = Utils.clear_xcorr(xcorr_train);
                                end

                                features_train = cat(2, mean_sd_train, xcorr_train, RPA_train, pos_reb_train, ROL_train, RPL_train, tLFP_train);
                                
                                T_train = cat(2, features_train, response);
                                assignin('base', 'T', T_train);

                                rng(7)
                                [trainedModel, valAcc] = trainClassifier(T_train, classnames, alg_type, alg);

                                validationAcc = cat(2, validationAcc, valAcc);

                                mean_sd_test = featureExtraction.mean_sd_grad_each_t(testData, 'n_submatrices', n_subm, 'issd', off_on{i_sd+1}, 'isgrad', off_on{i_grad+1});
                                mean_sd_test = Utils.clear_zeros_mean_sd_grad(mean_sd_test);
                                if i_xcorr == 0
                                    xcorr_test = [];
                                else
                                    xcorr_test = featureExtraction.xcorr2_all(testData);
                                    xcorr_test = Utils.clear_xcorr(xcorr_test);
                                end

                                T_test = cat(2, mean_sd_test, xcorr_test, RPA_test, pos_reb_test, ROL_test, RPL_test, tLFP_test);



                                % fit and predict
                                yfit_trainedModel = trainedModel.predictFcn(T_test); % fitting
                                Acc_trainedModel = Utils.accuracy(yfit_trainedModel, correct); % accuracy

                                testAcc= cat(2, testAcc, Acc_trainedModel);
                                fileID = fopen('all_features.txt','a');
                                fprintf(fileID,'n_subm %d, xcorr %d, issd %d, isgrad %d\n', n_subm, i_xcorr, i_sd, i_grad);
                                fprintf(fileID,'RPA %d, pos_reb %d, ROL %d, RPL %d, tLFP %d -- validationAcc %4.4f, testAcc %4.4f\n', i_RPA, i_pos_reb, i_ROL, i_RPL, i_tLFP, valAcc, Acc_trainedModel);
                                fclose(fileID);
                            end
                        end
                    end
                end   
            end
        end
        
        function confusion_matrix(yfit, correct, classnames, varargin)
            
            n_samples = length(yfit);
            targets = zeros(length(classnames), n_samples);
            outputs = targets;
            for i = 1:n_samples
                class = find(classnames==correct(i));
                targets(class, i) = 1;
                class = find(classnames==yfit(i));
                outputs(class, i) = 1;
            end
            plotconfusion(targets, outputs, varargin{:})
            xtickangle(0)
            if classnames(1) ~= 1
                xticklabels({'0.5V', '1.0V', '1.4V'})
                yticklabels({'0.5V', '1.0V', '1.4V'})
            end
            set(gca, 'fontsize', 12)
            %confusionmat
            % plotConfMat(targets, outputs)
%             if length(classnames) == 7
%                 for i = 1:n_samples
%                     targets(correct(i), i) = 1;
%                     outputs(yfit(i), i) = 1;
%                 end
%             else
%                 for i = 1:n_samples
%                    if correct(i) == 0.5
%                        targets(1, i) = 1;
%                    elseif correct(i) == 1.0
%                        targets(2,i) =1;
%                 end
%             end
            
        end

        function varargout = LDA_predict(MdlLinear, X)
            % Perform prediction using linear discriminant analysis
            % input:
            % - trained model
            % - feature matrix
            % output:
            % - predicted class
            % - score = posterior distribution for discriminant analysis in
            % Matlab, otherwise is likelihood.

            
            %cost = MdlLinear.Cost;
            %prior = MdlLinear.Prior;
            %sigma_k = MdlLinear.Sigma;
            %mu_k = MdlLinear.Mu;
            
            nSamples = size(X,1);
            nClass = length(MdlLinear.ClassNames);
            argmax = zeros(nSamples, 1);
            posterior = zeros(nSamples, nClass);
            
            if iscell(MdlLinear.ClassNames)
                predictedClass = cell(nSamples,1);
            else
                predictedClass = zeros(nSamples, 1);
            end
            
            mult_norm_const = (1/(sqrt(2*pi*det(MdlLinear.Sigma))));
            
            if ~isreal(mult_norm_const)
                mult_norm_const = -log(-imag(mult_norm_const));
            else
                mult_norm_const = log(mult_norm_const);
            end
            
            assignin('base', 'mult_norm_const', mult_norm_const);
            
            invSigma = pinv(MdlLinear.Sigma);
            for i = 1:nSamples
                post = zeros(nClass,1);
                post_cost = zeros(nClass,1);
                multivar_norm = zeros(nClass,1);
                norm_const = zeros(nClass,1);
                for k = 1:nClass
                    multivar_norm(k) = exp(mult_norm_const + ...
                        (-(X(i,:)-MdlLinear.Mu(k,:))*invSigma*(X(i,:)-MdlLinear.Mu(k,:))'/2));
                    
%               From book Bishop 1/nthrooth(2*pi, dimension of features space) is used
%                     multivar_norm(k) = (1/nthroot(2*pi, 26))*(1/(sqrt(det(sigma_k)))) * ...
%                         exp(-(test(i,:)-mu_k(k,:))*pinv(sigma_k)*(test(i,:)-mu_k(k,:))'/2);
                    norm_const(k) = multivar_norm(k)*MdlLinear.Prior(k);
                end
                
                assignin('base', 'multivar_norm', multivar_norm);
                assignin('base', 'norm_const', norm_const);
                for k = 1:nClass
                    post(k) = norm_const(k)/sum(norm_const);
                    post_cost(k) = sum(post(k)*MdlLinear.Cost(k,:));
                end
                assignin('base', 'post', post);
                
                [tmp_value(i), argmax(i)] = max(post_cost);
                posterior(i,:) = post;
                
                %[max_value(i) argmax(i)] = max(post_cost);
                %disp(posterior)
                % Note: on Matlab page explaining the prediction function it's
                % said that argmin is taken with minus sign in the argument of exp
                % in multivariate norm func,
                % instead here the argmax is taken with the minus sign. Argmin
                % works with plus sign in the argument of exp
            
                if iscell(MdlLinear.ClassNames)
                    predictedClass{i} = MdlLinear.ClassNames{argmax(i)};
                else
                    predictedClass(i) = MdlLinear.ClassNames(argmax(i));
                end
            end
            %assignin('caller', 'posterior', posterior);
            varargout{1} = predictedClass;
            varargout{2} = posterior;
            
            assignin('base', 'norm_const', norm_const)
        end
        function varargout = ESD_predict(trainedModel, T_test, varargin)
            % Prediction using ensemble subspace discriminant
            % Input:
            % - trained model
            % - feature matrix
            % - As third input if you insert 'vote', the prediction is based
            % on the voting of each cycle, i.e. the predicted class has the 
            % max number of votes, without this input the
            % prediction is made taking the argmax of mean score of all the
            % cycles. Default: argmax of mean score
            % output:
            % - predicted class
            % - score = mean posterior distribution for discriminant analysis in
            % Matlab, otherwise is likelihood.
            
            totCycle = trainedModel.ClassificationAlg.NumTrained;
            nClasses = length(trainedModel.ClassificationAlg.ClassNames);
            nSamples = size(T_test, 1);
            
            votes = zeros(nSamples, nClasses);
            predicted_classes = zeros(nSamples, 1);
            score = zeros(nSamples, nClasses, totCycle);
            for nCycle = 1:totCycle
                nFeatures = length(trainedModel.ClassificationAlg.PredictorNames);
                i = 1; test = [];
                for p = 1:nFeatures
                    if strcmp(trainedModel.ClassificationAlg.PredictorNames{p}, ...
                            trainedModel.ClassificationAlg.Trained{nCycle}.PredictorNames{i})
                        test = cat(2, test, T_test(:, p));
                        %disp(size(test))
                        i = i+1;
                    end
                    if i > length(trainedModel.ClassificationAlg.Trained{nCycle}.PredictorNames)
                        break
                    end
                end
                
                [predClass, posterior] = Utils.LDA_predict(trainedModel.ClassificationAlg.Trained{nCycle}, test);
                score(:,:, nCycle) = posterior;
                for i = 1:nSamples
                    if iscell(trainedModel.ClassificationAlg.ClassNames)
                        arg = find(predClass{i} == trainedModel.ClassificationAlg.ClassNames);
                    else
                        arg = find(predClass(i) == trainedModel.ClassificationAlg.ClassNames);
                    end
                    votes(i, arg) = votes(i, arg) + 1;
                end
            end
            %assignin('base', 'votes', votes)
            mean_score = mean(score,3);
            if isempty(varargin)
                % argmax mean score
                [M, argM] = max(mean_score, [],2);
                for i = 1:nSamples
                    predicted_classes(i) = trainedModel.ClassificationAlg.ClassNames(argM(i));
                end
            elseif strcmp(varargin{1}, 'vote')
                % voting system
                [max_votes, argmax_votes] = max(votes,[], 2);
                for i = 1:nSamples
                    predicted_classes(i) = trainedModel.ClassificationAlg.ClassNames(argmax_votes(i));
                end
            end
                
            %assignin('base','myScore', mean(score,3))
            
            varargout{1} = predicted_classes;
            varargout{2} = mean_score;
        end
        function compare_results(predicted_classes, yfit_trainedModel)
            compare = predicted_classes == yfit_trainedModel;
            fprintf('number of equal predictions: %d out of %d\n\n', sum(compare), length(yfit_trainedModel))
        end
        
        function pred_class = ESD_predict_C(trainedModel, T)
            TOT_CYCLES = 30;
            SUB_DIM = 49;
            N_CLASSES = 7;
            for cycle = 1:TOT_CYCLES
                Mu = trainedModel.ClassificationAlg.Trained{cycle}.Mu';
                semivector_Mu(cycle, :) = Mu(:)';
            end
            semivector_Mu = semivector_Mu';
            vector_Mu = semivector_Mu(:)';
            for cycle = 1:TOT_CYCLES
                pinv_sigma = pinv(trainedModel.ClassificationAlg.Trained{cycle}.Sigma(:,:));
                pinv_sigma = pinv_sigma';
                semivector_pinv_sigma(cycle, :) = pinv_sigma(:)';
            end
            semivector_pinv_sigma = semivector_pinv_sigma';
            vector_pinv_sigma = semivector_pinv_sigma(:)';
            for cycle = 1:TOT_CYCLES
                ln_k_det(cycle) = (1/(sqrt(2*pi*det(trainedModel.ClassificationAlg.Trained{cycle}.Sigma))));

                if ~isreal(ln_k_det(cycle))
                    ln_k_det(cycle) = -log(-imag(ln_k_det(cycle)));
                else
                    ln_k_det(cycle) = log(ln_k_det(cycle));
                end
            end
            [row, col] = find(trainedModel.ClassificationAlg.UsePredForLearner); % row
            subspace_idx = uint8(row)';
            
            acc_score = zeros(1, N_CLASSES);
            for cycle = 1:TOT_CYCLES
                for i_sub = 1:SUB_DIM
                    i_sub_idx = subspace_idx((cycle-1)*SUB_DIM + i_sub);
                    subX(i_sub) = T(i_sub_idx);
                end
                tmp_mu = vector_Mu((cycle-1)*N_CLASSES*SUB_DIM+1:cycle*N_CLASSES*SUB_DIM);
                tmp_pinv_sigma = vector_pinv_sigma((cycle-1)*SUB_DIM*SUB_DIM+1:cycle*SUB_DIM*SUB_DIM);
                lda_posterior = Utils.LDA_predict_C(subX, tmp_mu, tmp_pinv_sigma, ln_k_det(cycle));
                acc_score = acc_score + lda_posterior;
            end
            mean_score = acc_score / TOT_CYCLES;
            [value pred_class] = max(mean_score);
        end
        function posterior = LDA_predict_C(subX, mu, pinv_sigma, ln_k_det)
            N_CLASSES =7;
            SUB_DIM = 49;
            PRIOR = 1/N_CLASSES;
            pinv_sigma_mat = reshape(pinv_sigma, [SUB_DIM, SUB_DIM]);
            for k = 1:N_CLASSES
                tmp_mu_k = mu((k-1)*SUB_DIM+1:k*SUB_DIM);
                exponent = (subX-tmp_mu_k)*pinv_sigma_mat*(subX-tmp_mu_k)';
                multivar_norm = ln_k_det - exponent/2.0;
                multivar_norm_prior(k) = multivar_norm * PRIOR;
            end
            norm_const = sum(multivar_norm_prior);
            posterior = multivar_norm_prior; % / norm_const;
        end
        
        function save_model_for_C(trainedModel)
           % transfer the trained model to C header file through .csv
            nCycles = trainedModel.ClassificationAlg.NumTrained;
            % mu
            for cycle = 1:nCycles
                Mu = trainedModel.ClassificationAlg.Trained{cycle}.Mu';
                semivector_Mu(cycle, :) = Mu(:)';
            end
            semivector_Mu = semivector_Mu';
            vector_Mu = semivector_Mu(:)';
            csvwrite('Mu.dat',vector_Mu);

            % pinv_sigma
            for cycle = 1:nCycles
                pinv_sigma = pinv(trainedModel.ClassificationAlg.Trained{cycle}.Sigma(:,:));
                pinv_sigma = pinv_sigma';
                semivector_pinv_sigma(cycle, :) = pinv_sigma(:)';
            end
            semivector_pinv_sigma = semivector_pinv_sigma';
            vector_pinv_sigma = semivector_pinv_sigma(:)';
            assignin('base', 'vector_pinv_sigma', vector_pinv_sigma);
            csvwrite('pinvSigma.dat',vector_pinv_sigma);

            % ln_k_det
            for cycle = 1:nCycles
                ln_k_det(cycle) = (1/(sqrt(2*pi*det(trainedModel.ClassificationAlg.Trained{cycle}.Sigma))));

                if ~isreal(ln_k_det(cycle))
                    ln_k_det(cycle) = -log(-imag(ln_k_det(cycle)));
                else
                    ln_k_det(cycle) = log(ln_k_det(cycle));
                end
            end
            csvwrite('ln_k_det.dat', ln_k_det);
            
            % subspace_idx
            [row, col] = find(trainedModel.ClassificationAlg.UsePredForLearner); % row
            subspace_idx = uint8(row)'-1;
            csvwrite('subspace_idx.dat',subspace_idx);

        end
        
        function save_data_for_C(data)
            
%             % X
%             csvwrite('X.dat',X);
            if ndims(data) ~= 3
                data = squeeze(data);
            end
            
            % 3D data matrix to 1D array
            for t = 1:size(data, 3)
                matrix_t = data(:,:,t);
                matrix_t = matrix_t';
                semivector_data(t,:) = matrix_t(:)';
            end
            semivector_data = semivector_data';
            vector_data = semivector_data(:)';
            % assignin('base', 'vector_data', vector_data);
            csvwrite('vector_data.dat',vector_data);
            
        end
        
        function [M, f] = freq_spectrum(data, time)
            % input: data and time
            % output: magnitude and frequency
            
            data = squeeze(data);
            N = length(data);

            % Time specifications:
            dt = time(2) - time(1);	% seconds per sample
            Fs = 1/dt;	% samples per second

            % zero-padding
            % Z = zeros(30,1);
            % data = [data;Z];

            % Fourier Transform:
            D = fftshift(fft(data));
            M = abs(D);
            [v,i] = max(M);
            M = M(i:end);

            % Frequency specifications:
            dF = Fs/N;	% hertz
            f = 0:dF:Fs/2;	% hertz
        end
        
        function plot_new_features(data, time, layer, stim)
            % Plot the mean1 of each matrix and show the extracted new features
            figure
            [RP, RP_i] = featureExtraction.response_peak(data,time);
            [PR, PR_i] = featureExtraction.positive_rebound(data,time);
            minimum = featureExtraction.mean_sd_grad_each_t(data, 'n_submatrices', 1);
            plot(time, minimum, '.', 'markersize', 15)
            hold on
            t_stim_on = find(time == 0);
            if isempty(t_stim_on)
                disp('no stimulation onset found')
            else
                plot(time(t_stim_on), minimum(t_stim_on), 'or', 'markerfacecolor', 'r')
            end
            [RO, RO_i] = featureExtraction.response_onset(data,time);
            plot(time(RO_i), RO, 'om', 'linewidth', 1.5)
            plot(time(PR_i), PR, 'og', 'linewidth', 1.5)
            plot(time(RP_i), RP, 'ob', 'linewidth', 1.5)
            [tLFP, ctr_value, ctr_idx] = featureExtraction.time_norm_LFP(data, time);
            %plot(time(ctr_idx), ctr_value, 'ok')
            plot([time(ctr_idx), time(RO_i)], [ctr_value, RO], 'r', 'linewidth', 1)
            plot(time(RO_i:ctr_idx), minimum(RO_i:ctr_idx), 'r', 'linewidth', 1)
            xlabel('Time (ms)')
            ylabel('Response (mV)')
            % xlim([-20, 100])
            %title(['LFP - layer ' layer ' - Stim. ', stim])
            title(['LFP'])
            %title('Local field potential')
            set(gca,'fontsize',20)
            % figurepalette 
        end
        function out = clear_zeros_mean_sd_grad(T_train)
            % for stim classification when the cropping in time does not
            % exclude the time point of stimulation onset which is needed
            % to extract the new features, the old features,
            % i.e. mean, standar deviation and gradients, are 0, so we have
            % to remove these.
            
            z = zeros(size(T_train, 1), 1);
            c = 1;
            while c <= size(T_train, 2)
                if prod(T_train(:,c) == z)
                    T_train(:,c) = [];
                else
                    c = c+1;
                end
            end
            out = T_train;
        end
        function out = clear_xcorr(T_train)
            % for stim classification when the cropping in time does not
            % exclude the time point of stimulation onset which is needed
            % to extract the new features, the old features,
            % i.e. mean, standar deviation and gradients, are 0, so we have
            % to remove these.
            T_size = size(T_train, 1);
            z_x_xcorr = zeros(size(T_train, 1), 1)-12;
            z_y_xcorr = zeros(size(T_train, 1), 1)+12;
            c = 1;
            while c <= size(T_train, 2)
                if prod(T_train(:,c) == z_x_xcorr) || prod(T_train(:,c) == z_y_xcorr)
                    T_train(:,c) = [];
                else
                    c = c+1;
                end
            end
            out = T_train;
        end
        function [nc, sd] = find_max_acc(varargin)
            % [nc, sd] = find_max_acc(varargin)
            % give validation accuracy vector and test accuracy vector as
            % inputs, then the output is the indexes of the max validation
            % accuracy intersect with max test accuracy
            % if only one input vector is given, [nc, sd] are the indexes
            % of max accuracy in the given vector
            
            if length(varargin) == 2
                valAcc_param = varargin{1};
                testAcc_param = varargin{2};
                % search jointly
                max_val_idx = find(valAcc_param == max(max(valAcc_param)));
                max_test_idx = find(testAcc_param == max(max(testAcc_param)));
                abs_max_idx = intersect(max_val_idx, max_test_idx)

                ncol = size(valAcc_param, 1);
                sd = fix(abs_max_idx(1)/ncol)+1;
                nc = mod(abs_max_idx(1),ncol);

                valAcc_param(nc, sd)
                testAcc_param(nc, sd)
            
            elseif length(varargin) == 1
                valAcc_param = varargin{1};
                % search only in one
                [nc,sd] = find(valAcc_param == max(max(valAcc_param)));

                valAcc_param([nc, sd])
                for i = 1:length(nc)
                    disp('cv')
                    valAcc_param(nc(i), sd(i))
                    %disp('test')
                    %testAcc_param(nc(i), sd(i))
                end
            else
                error('too many inputs')
            end
            
        end
        
        function A = trapz(data)
            A = 0;
            for i = 1:length(data)-1
                A = A + (data(i)+data(i+1))/2.0;
            end

        end
    end
end