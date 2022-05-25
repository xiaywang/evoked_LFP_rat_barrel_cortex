%% plot lfp for ppt on 21.12.2017 for skype meeting with university of Padua and plots for the paper on layer classification

startingFolder = pwd;
foldername = uigetdir(startingFolder, 'Select the folder containing the data');
layers = ["0200", "0350", "0500", "0750", "1100", "1500", "1750"];

%layers = ["0180", "0360", "0540", "0720", "0900", "1440", "1620"];
%% extract one sample per layer
sample = zeros(7, 1, 16,16,69);
sample_mean1 = zeros(69, 7);
% figure()
for i = 1:length(layers)
    layer = layers(i);
    %subplot(7,1,i)
    %load(strcat(foldername, '/long-300ms-stim_amp1.4V-', layer, 'um-Group Name #1.mat'))
    load(strcat(foldername, '/short-5ms-stim_amp1.4V-', layer, 'um-Group Name #1.mat')); % short-5ms-stim_amp1.4V-0200um-Group Name #1
    %mean1 = (i-1)+featureExtraction.mean_sd_grad_each_t(data, 'n_submatrices', 1);
    % a(i) = plot(time, mean1, 'LineWidth', 1.5); m{i} = strcat('-', layer, 'um (L', num2str(i), ')');
    % ylim([-2, 2])
    sample(i,1,:,:,:) = data;
    mean1 = featureExtraction.mean_sd_grad_each_t(sample(i,1,:,:,:), 'n_submatrices', 1);
    %sample(:,:,:,i) = data;
    sample_mean1(:,i) = mean1;
end
% legend(a, m, 'Box', 'off', 'Position',[0.9 0.12 0.1 0.8]);
% xlabel('time (ms)')

%% plot 3d data

[y,z] = meshgrid(1:16,1:16);
x = ones(16, 16);
surf(x,y,z)
hold on
x = 2*ones(16, 16);
surf(x,y,z)

%% plot image data in 3d (showing time) and the average in time
layers_label = [200, 350, 500, 750, 1100, 1500, 1750];
t_start = 14;
t_end = 17;
bottom = min(min(min(min(sample(:, 1, :,:,t_start:t_end)))));
top  = max(max(max(max(sample(:, 1, :,:,t_start:t_end)))));
ylim_min = min(min(sample_mean1(t_start:t_end, :)));
ylim_max = max(max(sample_mean1(t_start:t_end, :)));

%<<<<<<< HEAD
%for i = 1:length(layers) %3:4
%    figure(3)
%    subplot(7,1,i)
%=======
for i = 7:7%1:length(layers)
    figure(1)
%>>>>>>> 6ce9d86d14b928e766835ce15e991a79585f300f
    %subplot(7,2,2*i-1)
    for t = t_start:t_end
        [y,z] = meshgrid(16:-1:1,16:-1:1);
        x = t*ones(16, 16);
        image = squeeze(sample(i,1,:,:,t));
        surf(x,y,z,image, 'edgecolor', 'none');
        hold on
        %contour3(x,y,z, [1,8,16], 'LineWidth',1, 'Color','g')
        %contour3(x,y,z, [1,16], 'LineWidth',1, 'Color','r')
    end
    caxis manual
    caxis([bottom top]);
    colormap(flipud(parula))
%<<<<<<< HEAD
%    colorbar
%    zlabel({[' Layer ', num2str(i), ' '], [num2str(layers_label(i)), 'um ']}, 'Rotation',0, 'HorizontalAlignment','right')
%=======
    c = colorbar
    c.Location = 'eastoutside';
%     c.Position = [0.9 0.1 0.3 0.7];
    %zlabel({[' Layer ', num2str(i), ' '], [num2str(layers_label(i)), 'um ']}, 'Rotation',0, 'HorizontalAlignment','right')
    zlabel('Recording site')
%>>>>>>> 6ce9d86d14b928e766835ce15e991a79585f300f
%     ax = gca;
%     ax.Visible = 'off';
%     if i < 7
%         ax = gca;
% %       ax.Visible = 'off';
% %         set(gca, 'zcolor', 'w')
%         set(gca, 'xticklabels', [])
%         set(gca, 'yticklabels', [])
% %         set(gca, 'zticklabels', [])
%         set(gca, 'xtick', [])
%         set(gca, 'ytick', [])
% %         set(gca, 'ztick', [])
%     end
    grid off
    set(gca, 'xticklabels', [])
%     xl = get(gca,'XLabel');
%     xlFontSize = get(xl,'FontSize');
%     disp(xlFontSize)
%     xAX = get(gca,'XAxis');
%     set(xAX,'FontSize', 10)
%     set(xl, 'FontSize', xlFontSize);
%     xlabel('Time')
    ylim([1,16])
    zlim([1,16])
    yticks([1 8 16])
    yticklabels({'16','8','1'})
    zticks([1 8 16])
    zticklabels({'16','8','1'})
%     xticks([14 15 16 17])
%     xticklabels({'11.2','12.8','14.4', '16'})
    if i == 7
        xlabel('Time')
%         ylabel({'Recording', 'site'})
        ylabel('Recording site')
    end
%     zlabel({'Recording site'})
    set(gca, 'fontsize', 22)
%     figurepalette
    
    figure(2)
    subplot(7,1,i)
    %subplot(7,2,2*i)
%<<<<<<< HEAD
%    plot(time, sample_mean1(:, i), 'LineWidth', 1)
%    %ylim([ylim_min-0.3, ylim_max+0.7])
%    ylim([-1.9, 0.9])
%=======
    plot(time, sample_mean1(:, i), 'LineWidth', 2)
    ylim([ylim_min-0.3, ylim_max+0.7])
%>>>>>>> 6ce9d86d14b928e766835ce15e991a79585f300f
    xlim([time(1)-1, time(end)])
    xlabel('Time (ms)')
    ylabel('Amplitude (mV)')
%     set(gca,'xtick',[],'ytick',[])
%     ax = gca;
%     ax.Visible = 'off';
    hold on
    for t = t_start:t_end
        plot(time(t), sample_mean1(t,i), 'ro', 'markersize', 6, 'markerfacecolor', 'r')
    end
    hold off
    if i < 7
      set(gca,'xticklabels',[]);
      set(gca, 'xlabel', [])
      set(gca, 'ylabel', [])
    end
    set(gca, 'fontsize', 22)
end
% set(gca,'xtickMode', 'auto')
% figurepalette


%% plot for features in paper

figure()
figurepalette
% subplot(1,2,1)
imshow(squeeze(sample(5, 1, :,:, 15)), [], 'Colormap', flipud(parula), 'InitialMagnification', 'fit')
c = colorbar;
c.FontSize = 40;
c.Ticks = -3:0.5:-0.5;
% xlabel('timepoint t')
% subplot(1,2,2)
figure()
figurepalette
imshow(squeeze(sample(5, 1, :,:, 16)), [], 'Colormap', flipud(parula), 'InitialMagnification', 'fit')
c = colorbar;
c.FontSize = 40;
c.Ticks = -3.5:0.5:-0.5;
% xlabel('timepoint t+1')

figure()
figurepalette
X = xcorr2(squeeze(sample(5, 1, 1+3:16-3,1+3:16-3, 15)), squeeze(sample(5, 1, :,:, 16)));
center = round(size(X)/2);
% subplot(2,2,[3,4]);
imshow(X, [], 'Colormap', flipud(parula), 'InitialMagnification', 'fit');
c = colorbar;
c.FontSize = 40;
hold on

[M_row, I_row] = max(X);
[M_col, I_col] = max(max(X));
if length(M_col) > 1
    disp('multiple max values')
    for m = 1:length(M_col)
        disp(find(X == M_col(m)))
    end
end
max_col = I_row(I_col);
max_row = I_col;
x = max_row-center(2);
y = center(1)-max_col;
plot(max_row, max_col, 'ro', 'Linewidth', 6)
quiver(center(2), center(1), x(end), -y(end), 'r','LineWidth', 8, 'MaxHeadSize', 1.5)
% title(strcat('xcorr t =', num2str(t), '-', num2str(t+1)));
%title('Cross-correlation map')
hold off

% featureExtraction.xcorr2_trial_layer(sample(:,:,17:18,4), 't_start', 1, 't_end', 2, 'plotting', 'on')

%%
x = 1:16; y = 1:16;
[X,Y] = meshgrid(x,y);

imshow(squeeze(sample(4, 1, :,:, 16)), [], 'Colormap', flipud(parula), 'InitialMagnification','fit')
hold on
%%
figure()
h = imagesc(squeeze(sample(4,1,:,:, 16)));
% set(gca, 'XTick', [1, 5, 10, 16]);
% set(gca, 'YTick', [1, 5, 10, 16]);
% xticklabels([1, 5, 10, 16])
% yticklabels([1, 5, 10, 16])
set(gca, 'XTick', [1, 8, 16]);
set(gca, 'YTick', [1, 8, 16]);
xticklabels([1, 8, 16])
yticklabels([1, 8, 16])
set(gca, 'FontSize', 28)
hold on;
for i = 1:16
   plot([.5,16.5],[i-.5,i-.5],'k-');
   plot([i-.5,i-.5],[.5,16.5],'k-');
end
plot([0.5 16.5], [16.5, 16.5], 'r-', 'linewidth', 10)
plot([0.5 16.5], [.5, .5], 'r-', 'linewidth', 10)
plot([16.5, 16.5], [0.5 16.5], 'r-', 'linewidth', 10)
plot([.5, .5], [0.5 16.5], 'r-', 'linewidth', 10)
% plot([.5, 16.5], [8.5 8.5], 'r-', 'linewidth', 10)
% plot([.5, 16.5], [5.5 5.5], 'r-', 'linewidth', 10)
% plot([.5, 16.5], [10.5 10.5], 'r-', 'linewidth', 10)
set(gcf, 'Position', [500, 300, 500, 500])
colormap(flipud(parula))

%% hyperparameters

% valAcc_tot = [];
% for nCycles = 5:100
%     rng(7);
%     [trainedModel, valAcc] = trainClassifier(T_train, 1:7, 'ensemble', 'nCycles', nCycles)
%     valAcc_tot = cat(1, valAcc_tot, valAcc);
% end

%close all
figure
subDim = 1:5:248;
valAcc_tot = [0.7746
    0.9018
    0.9576
    0.9732
    0.9844
    0.9888
    0.9933
    0.9933
    0.9955
    0.9911
    0.9888
    0.9911
    0.9933
    0.9933
    0.9955
    0.9955
    0.9933
    0.9911
    0.9933
    0.9866
    0.9866
    0.9933
    0.9866
    0.9844
    0.9911
    0.9866
    0.9888
    0.9754
    0.9866
    0.9911
    0.9799
    0.9844
    0.9866
    0.9821
    0.9844
    0.9799
    0.9821
    0.9732
    0.9777
    0.9710
    0.9732
    0.9754
    0.9732
    0.9710
    0.9688
    0.9754
    0.9643
    0.9710
    0.9754
    0.9643];
plot(subDim, 100*valAcc_tot, '.')
xlabel('Subspace Dimension')
ylabel('Validation Accuracy (%)')
xlim([0,250])
set(gca, 'fontsize', 16)
set(gcf, 'Position', [284.0000  448.5000  406.5000  319.5000])

figure()
nCycles = 5:100;
valAcc_tot = [0.9844
    0.9888
    0.9866
    0.9888
    0.9911
    0.9911
    0.9933
    0.9933
    0.9888
    0.9888
    0.9933
    0.9955
    0.9911
    0.9911
    0.9933
    0.9933
    0.9888
    0.9933
    0.9888
    0.9933
    0.9888
    0.9911
    0.9911
    0.9955
    0.9933
    0.9955
    0.9933
    0.9911
    0.9911
    0.9933
    0.9955
    0.9933
    0.9888
    0.9933
    0.9911
    0.9911
    0.9933
    0.9933
    0.9933
    0.9955
    0.9911
    0.9933
    0.9911
    0.9933
    0.9933
    0.9888
    0.9955
    0.9933
    0.9955
    0.9888
    0.9888
    0.9955
    0.9911
    0.9911
    0.9911
    0.9933
    0.9933
    0.9866
    0.9888
    0.9911
    0.9911
    0.9911
    0.9933
    0.9955
    0.9955
    0.9933
    0.9911
    0.9911
    0.9933
    0.9955
    0.9933
    0.9933
    0.9911
    0.9933
    0.9933
    0.9888
    0.9911
    0.9911
    0.9955
    0.9911
    0.9911
    0.9933
    0.9911
    0.9911
    0.9933
    0.9933
    0.9933
    0.9933
    0.9933
    0.9933
    0.9866
    0.9933
    0.9911
    0.9955
    0.9911
    0.9955];
plot(nCycles, 100*valAcc_tot, 'linewidth', 1.5)
xlabel('Number of Learners')
ylabel('Validation Accuracy (%)')
% yticks(95:0.5:100)
ylim([98, 100])
xticks(0:30:100)
set(gca, 'fontsize', 16)
set(gcf, 'Position', [284.0000  448.5000  406.5000  319.5000])

%% 

% calculate area using trapz
figure
x = 1:0.2:10;
y = x;
plot(x, y, '.')
(x(2) - x(1))*trapz(y)


%% plot new features

l = 3; % layer
s = 1; % sample
data = trainData(l, s, :,:,:);
Utils.plot_new_features(data, time, num2str(l), '1.4V');
figurepalette

%% explore new features (on new dataset)
% close all

[pr,i_pr] = featureExtraction.positive_rebound(trainData, time);
[rp,i_rp] = featureExtraction.response_peak(trainData, time);
[ro,i_ro] = featureExtraction.response_onset(trainData, time);

%%

%i_weird = find(i_rp == 25);
i_weird = 101;
stim = ceil(i_weird/nTrain);
sample = mod(i_weird, nTrain);
idx = (stim-1)*nTrain + sample;
mean1 = featureExtraction.mean_sd_grad_each_t(trainData(stim,sample,:,:,:), 'n_submatrices', 1);
figure
plot(time, mean1, '.')
hold on
plot(time(i_pr(idx)), pr(idx), 'ro')
plot(time(i_rp(idx)), rp(idx), 'bo')
plot(time(i_ro(idx)), ro(idx), 'go')

%%
figure
subplot(121)
histogram(pr)
subplot(122)
histogram(i_pr)
title('positive rebound')

figure
subplot(121)
histogram(rp)
subplot(122)
histogram(i_rp)
title('response peak')

figure
subplot(121)
histogram(ro)
subplot(122)
histogram(i_ro)
title('response onset')

%%

i_weird = 101;
stim = ceil(i_weird/nTrain);
sample = mod(i_weird, nTrain);
data = trainData(stim,sample,:,:,:);
Utils.plot_new_features(data, time, num2str(l), [num2str(stim) 'V']);
figurepalette



%% plot_xcorr_feature_analysis

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

%% training set

% preprocessing
[trainData, time] = Preprocessing.cropData(trainData, time, 8, 32);
sigma = 0.6;
trainData = Preprocessing.gaussian_smooth(trainData, sigma);

%%
xcorr_train = featureExtraction.xcorr2_all(trainData);

x_corr = xcorr_train(:,1:24);
y_corr = xcorr_train(:, 25:48);

x_corr_l = zeros(7, 64*24);
y_corr_l = zeros(7, 64*24);
for l = 1:7
    tmp = x_corr((l-1)*64+1:l*64, :);
    x_corr_l(l,:) = reshape(tmp, [1,64*24]);
    tmp = y_corr((l-1)*64+1:l*64, :);
    y_corr_l(l,:) = reshape(tmp, [1,64*24]);
end

%% modulus

modulus = sqrt(x_corr_l.^2+y_corr_l.^2);
modulus = modulus';
angle = radtodeg(atan(y_corr_l./x_corr_l));
angle = angle';
modulus(find(modulus==0)) = NaN;

f1 = figure('name', 'Modulus', 'NumberTitle', 'off');
boxplot(modulus, 'Whisker',1.5)
set(findobj(gca,'type','line'),'linew',1)
xlabel('Layers')
ylabel('Modulus')
title('Modulus - data set 1.4V')
set(gca, 'fontsize', 18)
f2 = figure('name', 'Angle', 'NumberTitle', 'off');
boxplot(angle, 'Whisker',1.5)
set(findobj(gca,'type','line'),'linew',1)
xlabel('Layers')
ylabel('Angle (deg)')
title('Angle - data set 1.4V')
set(gca, 'fontsize', 18)

mean_modulus = nanmean(modulus);
sd_modulus = nanstd(modulus);
mean_angle = nanmean(angle);
sd_angle = nanstd(angle);

%%
figure
histogram(modulus(:,2))
set(gca, 'fontsize', 14)
title('Modulus of layer 2')
xlabel('Modulus')
ylabel('Counts')

figure
histogram(angle(:,2))
set(gca, 'fontsize', 14)
title('Angle of layer 2')
xlabel('Angle')
ylabel('Counts')

%% RPA histograms

% merged data
startingFolder = pwd;
foldername = uigetdir(startingFolder, 'Select the folder containing all the different data folders');
% layers = ["0240", "0360", "0480", "0720", "0960", "1440", "1680"];
% layers = ["0500", "0750", "1750"];
layers = ["0200", "0350", "0500", "0750", "1100", "1500", "1750"];
%layers = ["0500"];
stimulations = [0.6; 0.8; 1.0; 1.2; 1.4];
%stimulations = [0.5; 1.0; 1.4];
[trainData, testData, time] = Preprocessing.merge(foldername, layers(3), 0.8, 0.2);

nTrain = size(trainData, 2); % number of samples for each class
nTest = size(testData, 2); % number of samples for each class

sigma = 0.6;
[trainData, time] = Preprocessing.cropData(trainData, time, 7, 45); %50
trainData = Preprocessing.gaussian_smooth(trainData, sigma); % HERE

[testData, time] = Preprocessing.cropData(testData, time, 7, 45); %50
testData = Preprocessing.gaussian_smooth(testData, sigma); % HERE

%%
NS = length(stimulations);

RPA_train = featureExtraction.response_peak_amplitude(trainData,time);
RPA_test = featureExtraction.response_peak_amplitude(testData,time);
tLFP_train = featureExtraction.time_norm_LFP(trainData, time);
tLFP_test = featureExtraction.time_norm_LFP(testData, time);
PR_train = featureExtraction.positive_rebound(trainData, time);
PR_test = featureExtraction.positive_rebound(testData, time);

RPA = zeros(NS, 80);
tLFP = zeros(NS, 80);
PR = zeros(NS, 80);
for i = 1:NS
    
    RPA_train_tmp = RPA_train((i-1)*(nTrain)+1:(i)*(nTrain));
    RPA_test_tmp = RPA_test((i-1)*(nTest)+1:(i)*(nTest));
    RPA(i, :) = [RPA_train_tmp; RPA_test_tmp];
    tLFP_train_tmp = tLFP_train((i-1)*(nTrain)+1:(i)*(nTrain));
    tLFP_test_tmp = tLFP_test((i-1)*(nTest)+1:(i)*(nTest));
    tLFP(i, :) = [tLFP_train_tmp; tLFP_test_tmp];
    PR_train_tmp = RPL_train((i-1)*(nTrain)+1:(i)*(nTrain));
    PR_test_tmp = RPL_test((i-1)*(nTest)+1:(i)*(nTest));
    PR(i, :) = [PR_train_tmp; PR_test_tmp];
    
end


%% histogram of RPA

figure
for i = 1:NS
    
    subplot(NS,1,i)
    histogram(RPA(i,:), 'binwidth', 1/40)
    xlim([0,0.5])
    xlabel('RPA(mV)')
    ylim([0, 35])
    ylabel('Counts')
    title(strcat(num2str(stimulations(i)), 'V'))
    set(gca, 'fontsize', 12)
    figurepalette
end

mean(RPA, 2)
std(RPA,[], 2)
median(RPA, 2)

%%
F = RPA;
figure
for i = 1:NS
    for j = 1:NS
        
        subplot(NS, NS, (i-1)*NS+j)
        
        scatter(F(i,:), F(j,:), 10, 'filled')
%         xlim([0,1])
%         ylim([0,1])
%         xticks(0:0.2:1.4)
%         yticks(0:0.2:1.4)
        pbaspect([1 1 1])
        if mod((i-1)*NS+j, NS) == 1
            ylabel(strcat(num2str(stimulations(i)), 'V'), 'fontweight','bold')
%             set(get(gca,'ylabel'),'rotation',0)
        end
        if (i-1)*NS+j < NS+1
            title(strcat(num2str(stimulations(j)), 'V'))
%             set(get(gca,'ylabel'),'rotation',0)
        end
    end
end

%% t-test 

% qq-plot
figure
for i = 1:NS
   subplot(3,2,i)
   qqplot(RPA(i,:))
   set(gca, 'fontsize', 12)
end

%%
decision = zeros(NS);
p_value = zeros(NS);
for i = 1:NS
    for j = 1:NS
        [decision(i,j), p_value(i,j)] = ttest2(RPA(i,:), RPA(j,:));
    end
end

%% SCATTER PLOTS

% normalize the features among all stimuli per each feature
% for i = 1:NS
%     for s = 1:length(RPA)
%         RPA(i,s) = (RPA(i,s) - min(RPA(i,:)))/(max(RPA(i,:)) - min(RPA(i,:)));
%         tLFP(i,s) = (tLFP(i,s) - min(tLFP(i,:)))/(max(tLFP(i,:)) - min(tLFP(i,:)));
%         RPL(i,s) = (RPL(i,s) - min(RPL(i,:)))/(max(RPL(i,:)) - min(RPL(i,:)));
%     end
% end

RPA(:,:) = (RPA(:,:) - min(min(RPA(:,:))))/(max(max(RPA(:,:))) - min(min(RPA(:,:))));
tLFP(:,:) = (tLFP(:,:) - min(min(tLFP(:,:))))/(max(max(tLFP(:,:))) - min(min(tLFP(:,:))));
PR(:,:) = (PR(:,:) - min(min(PR(:,:))))/(max(max(PR(:,:))) - min(min(PR(:,:))));

myColors = zeros(length(RPA) + length(tLFP) + length(PR), 3);
for i = 1:80
    myColors(i, :) = [1, 0, 0]; % red
    myColors(80+i, :) = [0, 1, 0]; % green
    myColors(160+i, :) = [0, 0, 1]; % blue
end
figure
for i = 1:NS
    for j = 1:NS
        
        subplot(NS, NS, (i-1)*NS+j)
        y = [RPA(i,:), tLFP(i,:), PR(i,:)];
        x = [RPA(j,:), tLFP(j,:), PR(j,:)];
%         x = [tLFP(i,:), RPA(i,:)];
%         y = [tLFP(j,:), RPA(j,:)];
        scatter(x, y, 8, myColors, 'filled')
        set(gca, 'fontsize', 14)
%         xlim([0,1])
%         ylim([0,1])
%         xticks(0:0.2:1.4)
%         yticks(0:0.2:1.4)
        pbaspect([1 1 1])
        if mod((i-1)*NS+j, NS) == 1
            ylabel(strcat(num2str(stimulations(i)), 'V'), 'fontweight','bold')
%             set(get(gca,'ylabel'),'rotation',0)
        end
        if (i-1)*NS+j < NS+1
            title(strcat(num2str(stimulations(j)), 'V'))
%             set(get(gca,'ylabel'),'rotation',0)
        end
    end
end

%% SPECTRUM

data06 = trainData(1,1,:);
data06 = squeeze(data06);
% data06 = featureExtraction.mean_sd_grad_each_t(trainData(1:16,:,:), 'n_submatrices', 1);
N = length(data06);

% Time specifications:
dt = (time_full(2) - time_full(1)) * 10^(-3);	% seconds per sample
Fs = 1/dt;	% samples per second


% zero-padding
% Z = zeros(30,1);
% data06 = [data06;Z];

% Fourier Transform:
D = fftshift(fft(data06));
M = abs(D);
[v,i] = max(M);
M = M(i:end);

% Frequency specifications:
dF = Fs/N;	% hertz
f = 0:dF:Fs/2;	% hertz

% Plot the spectrum:
figure;
plot(f,M);
xlabel('Frequency [Hz]');
title('Magnitude Response');

%% low pass FIR filter


   
%%

data06 = trainData(1,1,:);
data06 = squeeze(data06);

D = fft(data06);
D = fftshift(D);

N = length(data06);
f = [-N/2:N/2-1]/N;

figure(1)
subplot(3,1,1)
plot(time_full, data06, '.')
subplot(3,1,2)
plot(f,abs(D))
% subplot(3,1,3)
% plot()

%%

data08 = trainData(17,1,:);
data08 = squeeze(data08);

D = fft(data08);
D = fftshift(D);

N = length(data08);
f = [-N/2:N/2-1]/N;

figure(2)
subplot(3,1,1)
plot(time_full, data08, '.-')
subplot(3,1,2)
plot(f,abs(D), '.-')
% subplot(3,1,3)
% plot()


%%

for s = 1:1024:5120
    
        data = trainData(1,1,:);

        [M,f] = Utils.freq_spectrum(data, time_full);

        % Plot the spectrum:
        figure;
        plot(f,M);
        xlabel('Frequency [Hz]');
        title('Magnitude Response');
    
end
