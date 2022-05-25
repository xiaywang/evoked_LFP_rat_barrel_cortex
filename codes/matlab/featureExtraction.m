classdef featureExtraction
    methods (Static)
        function out = mean_sd_grad_each_t(inputData, varargin)
            % extract features concatenating them by columns and giving them as output
            % By default it computes the mean of half matrix (no std, no gradient)
            % It's possible to compute mean and sd of whole matrix; mean
            % std and gradient (mean of top half matrix minus mean of
            % bottom half matrix) of half matrices; mean and std of 1/3 and
            % 1/4 matrices.
            % Use the following to set the desired features:
            % - 'n_submatrices', to divide the 16x16 matrix into submatrices to compute the features,
            % 1 means no division (i.e. whole matrix), 2 means top and
            % bottom half matrices, 3 means top, central and bottom
            % matrices, 4 means top left, top right, bottom left and bottom
            % right (for 4 you can set overlap 'on' if you want to have top half, bottom half, left half and right half)
            % - 'issd', set it to 'on' if you want to compute std (default 'off')
            % - 'isgrad', set it to 'on' if you want to compute the gradient
            % of means (default 'off')
            
            n_submatrices = 2; issd = 'off'; isgrad = 'off'; overlap = 'off';
            
            for i = 1:length(varargin)
                if strcmp(varargin{i}, 'n_submatrices')
                    n_submatrices = varargin{i+1};
                    if n_submatrices > 4 || n_submatrices < 1
                        error('Incorrect n_submatrices, possible values: 1, 2, 3, 4')
                    end
                elseif strcmp(varargin{i}, 'issd')
                    issd = varargin{i+1};
                elseif strcmp(varargin{i}, 'isgrad')
                    isgrad = varargin{i+1};
                elseif strcmp(varargin{i}, 'overlap')
                    overlap = varargin{i+1};
                end
            end
            
            out1 = []; out2 = []; out3 = []; out4 = [];
            if n_submatrices == 1
                out1 = featureExtraction.mean1_sd1_each_t(inputData, 'issd', issd);
            elseif n_submatrices == 2
                out2 = featureExtraction.mean2_sd2_each_t(inputData, 'issd', issd, 'isgrad', isgrad);
            elseif n_submatrices == 3
                out3 = featureExtraction.mean3_sd3_each_t(inputData, 'issd', issd, 'isgrad', isgrad);
            elseif n_submatrices == 4
                out4 = featureExtraction.mean4_sd4_each_t(inputData, 'issd', issd, 'isgrad', isgrad, 'overlap', overlap);
            else
                error('Invalid n_submatrices')
            end
            out = cat(2, out1, out2, out3, out4);
        end
        
        function out = mean1_sd1_each_t(inputData, varargin)
            % Takes the data and gives the mean of 16-by-16 matrix. If the
            % name 'issd' is given as additional input followed by 'on',
            % then also the stadard deviation is given per each time point.
            % Default is 'off'.
            % input:
            %   5D data: (layer, n_samples, x, y, time)
            %   varargin: 'issd', 'off' --> it doesn't give the standard
            %   deviation. 'issd', 'on' --> it performs mean and standard
            %   deviation
            % output:
            %   mean of 16-by-16 matrix per each time point or mean and standard
            %   deviation of 16-by-16 matrix per each time point
            
            issd = 'off';
            for i = 1:length(varargin)
                if strcmp(varargin{i}, 'issd')
                    issd = varargin{i+1};
                end
            end
            
            out = [];
            for t = 1:size(inputData, 5)
                mean_1 = []; sd_1 = [];
                for l = 1:size(inputData, 1) % l = layer
                    for i = 1:size(inputData, 2)
                        mean_1 = cat(1, mean_1, mean(mean(inputData(l,i,:,:,t))));
                        if strcmp(issd, 'on')
                            sd_1 = cat(1, sd_1, std2(inputData(l,i,:,:,t)));
                        end
                    end
                end
                if strcmp(issd, 'on')
                    out = cat(2, out, mean_1, sd_1);
                else
                    out = cat(2, out, mean_1);
                end
            end
        end
        function out = mean2_sd2_each_t(inputData, varargin)
            
            issd = 'off'; isgrad = 'off';
            
            for i = 1:length(varargin)
                if strcmp(varargin{i}, 'issd')
                    issd = varargin{i+1};
                elseif strcmp(varargin{i}, 'isgrad')
                    isgrad = varargin{i+1};
                end
            end
            
            out = [];
            for t = 1:size(inputData, 5)
                mean_up = []; mean_bottom = []; sd_up = []; sd_bottom = [];
                for l = 1:size(inputData, 1)
                    for i = 1:size(inputData, 2)
                        bottom = size(inputData, 3);
                        center = bottom/2;
                        mean_up = cat(1, mean_up, mean(mean(inputData(l,i,1:center,:,t))));
                        mean_bottom = cat(1, mean_bottom, mean(mean(inputData(l,i,center+1:bottom,:,t))));
                        if strcmp(issd, 'on')
                            sd_up = cat(1, sd_up, std2(inputData(l,i,1:center,:,t)));
                            sd_bottom = cat(1, sd_bottom, std2(inputData(l,i,center+1:bottom,:,t)));
                        end
                    end
                end
                if strcmp(issd, 'on') && strcmp(isgrad, 'on')
                    out = cat(2, out, mean_up, mean_bottom, sd_up, sd_bottom, mean_up-mean_bottom);
                elseif strcmp(issd, 'off') && strcmp(isgrad, 'on')
                    out = cat(2, out, mean_up, mean_bottom, mean_up-mean_bottom);
                elseif strcmp(issd, 'on') && strcmp(isgrad, 'off')
                    out = cat(2, out, mean_up, mean_bottom, sd_up, sd_bottom);
                elseif strcmp(issd, 'off') && strcmp(isgrad, 'off')
                    out = cat(2, out, mean_up, mean_bottom);
                end
            end
        end
        function out = mean3_sd3_each_t(inputData, varargin)
            
            issd = 'off'; isgrad = 'off';
            
            for i = 1:length(varargin)
                if strcmp(varargin{i}, 'issd')
                    issd = varargin{i+1};
                elseif strcmp(varargin{i}, 'isgrad')
                    isgrad = varargin{i+1};
                end
            end
            
            out = [];
            for t = 1:size(inputData, 5)
                mean_up = []; mean_central = []; mean_bottom = [];
                sd_up = []; sd_central = []; sd_bottom = [];
                for l = 1:size(inputData, 1)
                    for i = 1:size(inputData, 2)
                        bottom = size(inputData, 3);
                        center1 = 5; center2 = 10;
                        mean_up = cat(1, mean_up, mean(mean(inputData(l,i,1:center1,:,t))));
                        mean_central = cat(1, mean_central, mean(mean(inputData(l,i,center1+1:center2,:,t))));
                        mean_bottom = cat(1, mean_bottom, mean(mean(inputData(l,i,center2+1:bottom,:,t))));
                        if strcmp(issd, 'on')
                            sd_up = cat(1, sd_up, std2(inputData(l,i,1:center1,:,t)));
                            sd_central = cat(1, sd_central, std2(inputData(l,i,center1+1:center2,:,t)));
                            sd_bottom = cat(1, sd_bottom, std2(inputData(l,i,center2+1:bottom,:,t)));
                        end
                    end
                end
                if strcmp(issd, 'on') && strcmp(isgrad, 'on')
                    out = cat(2, out, mean_up, mean_central, mean_bottom, sd_up, sd_central, sd_bottom, mean_up-mean_central, mean_central-mean_bottom);
                elseif strcmp(issd, 'off') && strcmp(isgrad, 'on')
                    out = cat(2, out, mean_up, mean_central, mean_bottom, mean_up-mean_central, mean_central-mean_bottom);
                elseif strcmp(issd, 'on') && strcmp(isgrad, 'off')
                    out = cat(2, out, mean_up, mean_central, mean_bottom, sd_up, sd_central, sd_bottom);
                elseif strcmp(issd, 'off') && strcmp(isgrad, 'off')
                    out = cat(2, out, mean_up, mean_central, mean_bottom);
                end
            end
        end
        function out = mean4_sd4_each_t(inputData, varargin)
            
            issd = 'off'; isgrad = 'off'; overlap = 'off';
            
            for i = 1:length(varargin)
                if strcmp(varargin{i}, 'issd')
                    issd = varargin{i+1};
                elseif strcmp(varargin{i}, 'isgrad')
                    isgrad = varargin{i+1};
                elseif strcmp(varargin{i}, 'overlap')
                    overlap = varargin{i+1};
                end
            end
            
            out = [];
            for t = 1:size(inputData, 5)
                mean_1 = []; mean_2 = []; mean_3 = []; mean_4 = [];
                sd_1 = []; sd_2 = []; sd_3 = []; sd_4 = [];
                for l = 1:size(inputData, 1)
                    for i = 1:size(inputData, 2)
                        bottom = size(inputData, 3); center = bottom/2;
                        if strcmp(overlap, 'off')
                            mean_1 = cat(1, mean_1, mean(mean(inputData(l,i,1:center,1:8,t))));
                            mean_2 = cat(1, mean_2, mean(mean(inputData(l,i,1:center,9:16,t))));
                            mean_3 = cat(1, mean_3, mean(mean(inputData(l,i,center+1:bottom,1:8,t))));
                            mean_4 = cat(1, mean_4, mean(mean(inputData(l,i,center+1:bottom,9:16,t))));
                            if issd > 0
                                sd_1 = cat(1, sd_1, std2(inputData(l,i,1:center,1:8,t)));
                                sd_2 = cat(1, sd_2, std2(inputData(l,i,1:center,9:16,t)));
                                sd_3 = cat(1, sd_3, std2(inputData(l,i,center+1:bottom,1:8,t)));
                                sd_4 = cat(1, sd_4, std2(inputData(l,i,center+1:bottom,9:16,t)));
                            end
                        elseif strcmp(overlap, 'on')
                            mean_1 = cat(1, mean_1, mean(mean(inputData(l,i,1:center,:,t))));
                            mean_2 = cat(1, mean_2, mean(mean(inputData(l,i,center+1:bottom,:,t))));
                            mean_3 = cat(1, mean_3, mean(mean(inputData(l,i,1:bottom,1:8,t))));
                            mean_4 = cat(1, mean_4, mean(mean(inputData(l,i,1:bottom,9:16,t))));
                            if strcmp(issd, 'on')
                                sd_1 = cat(1, sd_1, std2(inputData(l,i,1:center,:,t)));
                                sd_2 = cat(1, sd_2, std2(inputData(l,i,center+1:bottom,:,t)));
                                sd_3 = cat(1, sd_3, std2(inputData(l,i,1:bottom,1:8,t)));
                                sd_4 = cat(1, sd_4, std2(inputData(l,i,1:bottom,9:16,t)));
                            end
                        end
                    end
                end
                if strcmp(issd, 'on') && strcmp(isgrad, 'on')
                    out = cat(2, out, mean_1, mean_2, mean_3, mean_4, sd_1, sd_2, sd_3, sd_4, mean_1-mean_2, mean_1-mean_3, mean_2-mean_4, mean_3-mean_4);
                elseif strcmp(issd, 'off') && strcmp(isgrad, 'on')
                    out = cat(2, out, mean_1, mean_2, mean_3, mean_4, mean_1-mean_2, mean_1-mean_3, mean_2-mean_4, mean_3-mean_4);
                elseif strcmp(issd, 'on') && strcmp(isgrad, 'off')
                    out = cat(2, out, mean_1, mean_2, mean_3, mean_4, sd_1, sd_2, sd_3, sd_4);
                elseif strcmp(issd, 'off') && strcmp(isgrad, 'off')
                    out = cat(2, out, mean_1, mean_2, mean_3, mean_4);
                end
            end
        end
        
        function out = xcorr2_all(data, varargin)
            if isempty(varargin)
                w = 10;
            else
                w = varargin{1};
            end
            
            x_corr = []; y_corr = [];
            crop = (16-w)/2;
            for l = 1:size(data, 1)
                for i = 1:size(data, 2)
                    x = []; y = [];
                    for t = 1:size(data, 5)-1
                        X = xcorr2(squeeze(data(l, i, 1+crop:16-crop, 1+crop:16-crop, t)), squeeze(data(l, i, :, :, t+1)));
    %                     X = xcorr2(data(top:bottom, :, t+1), data(top+crop:bottom-crop, 1+crop:16-crop, t));
                        center = round(size(X)/2);
                        [M_row, I_row] = max(X);
                        [M_col, I_col] = max(max(X));
                        if length(M_col) > 1 % in case of multiple values
                            disp('multiple max values')
                            for m = 1:length(M_col)
                                disp(find(X == M_col(m)))
                            end
                        end
                        max_col = I_row(I_col);
                        max_row = I_col;
                        x = cat(2, x, max_row-center(2));
                        y = cat(2, y, center(1)-max_col);
                    end
                    x_corr = cat(1, x_corr, x);
                    y_corr = cat(1, y_corr, y);
                end
            end
            out = cat(2, x_corr, y_corr);
        end
        
        function [x, y] = xcorr2_trial_layer(data, varargin)
            % cross-correlation between subsequent frames using the first
            % cropped frame (w-by-w) as kernel to correlate with the
            % subsequent frame.
            % 'trial' = number of trial (default 1)
            % 'layer' = number of layer (default 1)
            % 't_start' and 't_end': cropping in time (default 8 and 30, respectively)
            % 'w' = window size - only even numbers (default 10)
            % 'plotting': on or off to show the cross-correlation and the frames
            % (default 'off')
            
            % output are the x and y coordinate of the vector from the
            % central point of the cross-correlation matrix to the higherst
            % value of cross-correlation
            
            
            trial = 1; layer = 1; t_start = 8; t_end = 30; w = 10; plotting = 'off';
            
            for i = 1:length(varargin)
                if strcmp(varargin{i}, 'trial')
                    trial = varargin{i+1};
                elseif strcmp(varargin{i}, 'layer')
                    layer = varargin{i+1};
                elseif strcmp(varargin{i}, 't_start')
                    t_start = varargin{i+1};
                elseif strcmp(varargin{i}, 't_end')
                    t_end = varargin{i+1};
                elseif strcmp(varargin{i}, 'w')
                    w = varargin{i+1};
                    if mod(w, 2)
                        error('Window size w is not even')
                        return
                    end
                elseif strcmp(varargin{i}, 'plotting')
                    plotting = varargin{i+1};
                end
            end
            
            x = []; y = [];
            nTrial = size(data, 1)/(16*7);
            up = (16*nTrial*(layer-1)) + (16*trial - 15);
            down = up + 15;
            t = t_start;
            crop = (16-w)/2;
            n_plots = t_end-t_start;
            if strcmp(plotting, 'on')
                nrows = floor(sqrt(n_plots));
                ncols = ceil(n_plots/nrows);
                figure
            end
            for p = 1:n_plots
                X = xcorr2(data(up+crop:down-crop, 1+crop:16-crop, t), data(up:down, :, t+1));
%                 X = xcorr2(data(up:down, :, t+1), data(up+crop:down-crop, 1+crop:16-crop, t));
                center = round(size(X)/2);
                if strcmp(plotting, 'on')
                    subplot(nrows, ncols, p)
                    imshow(X, [], 'Colormap', flipud(parula), 'InitialMagnification', 'fit');
                    colorbar
                    hold on
                end
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
                x = cat(1, x, max_row-center(2));
                y = cat(1, y, center(1)-max_col);
                if strcmp(plotting, 'on')
                    plot(max_row, max_col, 'ro')
                    quiver(center(2), center(1), x(end), -y(end), 'r','LineWidth', 1.5, 'MaxHeadSize', 1.5)
                    title(strcat('xcorr t =', num2str(t), '-', num2str(t+1)));
                    hold off
                end
                t = t+1;
            end
            if strcmp(plotting, 'on')
                h.NextPlot = 'add';
                a = axes; 
                %// Set the title and get the handle to it
                ht = title(strcat('trial ', num2str(trial), ' - layer ', num2str(layer)));
                %// Turn the visibility of the axes off
                a.Visible = 'off';
                %// Turn the visibility of the title on
                ht.Visible = 'on';
                figure
                if n_plots+1 > nrows*ncols
                    ncols = ncols+1;
                end
                t = t_start;
                for p = 1:n_plots+1
                    subplot(nrows, ncols, p)
                    imshow(data(up:down, :, t), [], 'Colormap', flipud(parula), 'InitialMagnification', 'fit');
                    colorbar
                    title(strcat('t=', num2str(t)));
                    % xlabel(num2str(mean(mean(data(up:down, :, t)))))
                    t = t+1;
                end
                h.NextPlot = 'add';
                a = axes; 
                %// Set the title and get the handle to it
                ht = title(strcat('trial ', num2str(trial), ' - layer ', num2str(layer)));
                %// Turn the visibility of the axes off
                a.Visible = 'off';
                %// Turn the visibility of the title on
                ht.Visible = 'on';
            end
        end
        function [U, V] = OpticalFlow_trial_layer(data, varargin)
            % Optical flow between subsequent frames using a kernel window
            % defined by ww
            % 'trial' = number of trial (default 1)
            % 'layer' = number of layer (default 1)
            % 't_start' and 't_end': cropping in time (default 15 and 16, respectively)
            % 'window_size' = ww (default 10)
            % 'plotting': on or off to show the cross-correlation and the frames
            % (default 'off')
            
            % output are the U and V coordinate of the optical flow vectors
            
            
            trial = 1; layer = 1; t_start = 15; t_end = 16; ww = 10; plotting = 'off';
            
            for i = 1:length(varargin)
                if strcmp(varargin{i}, 'trial')
                    trial = varargin{i+1};
                elseif strcmp(varargin{i}, 'layer')
                    layer = varargin{i+1};
                elseif strcmp(varargin{i}, 't_start')
                    t_start = varargin{i+1};
                elseif strcmp(varargin{i}, 't_end')
                    t_end = varargin{i+1};
                elseif strcmp(varargin{i}, 'window_size')
                    ww = varargin{i+1};
                    if mod(ww, 2)
                        error('Window size w is not even')
                        return
                    end
                elseif strcmp(varargin{i}, 'plotting')
                    plotting = varargin{i+1};
                end
            end
            
            U = []; V = [];
            
            nTrial = size(data, 1)/(16*7);
            up = (16*nTrial*layer) + (16*trial - 15);
            down = up + 15;
            t = t_start;
            n_plots = t_end-t_start;
            if strcmp(plotting, 'on')
                nrows = floor(sqrt(n_plots));
                ncols = ceil(n_plots/nrows);
                figure
                fr1 = data(up:down,:,t);
                fr2 = data(up:down,:,t+1);
                subplot 121
                imshow(fr1, [], 'Colormap', flipud(parula), 'InitialMagnification','fit');
                subplot 122
                imshow(fr2, [], 'Colormap', flipud(parula), 'InitialMagnification','fit');

                figure
            end
            for p = 1:n_plots
                fr1 = imresize(data(up:down,:,t), [64, 64]);
                fr2 = imresize(data(up:down,:,t+1), [64, 64]);
                [X, Y, u, v] = featureExtraction.OF(fr1, fr2, ww);
                U = cat(1, U, u);
                V = cat(1, V, v);
                if strcmp(plotting, 'on')
                    subplot(nrows, ncols, p)
                    imshow(fr2, [], 'Colormap', flipud(parula), 'InitialMagnification','fit');
                    hold on;
                    % draw the velocity vectors
                    quiver(X, Y, u, v, 'w', 'LineWidth', 0.8, 'MaxHeadSize', 1)
                    xlabel(strcat('Optical Flow t =', num2str(t), '-', num2str(t+1)));
                    hold off
                end
                t = t+1;
            end
            if strcmp(plotting, 'on')
                h.NextPlot = 'add';
                a = axes; 
                %// Set the title and get the handle to it
                ht = title(strcat('trial ', num2str(trial), ' - layer ', num2str(layer)));
                %// Turn the visibility of the axes off
                a.Visible = 'off';
                %// Turn the visibility of the title on
                ht.Visible = 'on';
                figure
                if n_plots+1 > nrows*ncols
                    ncols = ncols+1;
                end
                t = t_start;
                for p = 1:n_plots+1
                    subplot(nrows, ncols, p)
                    fr = imresize(data(up:down,:,t), [64, 64]);
                    imshow(fr, [], 'Colormap', flipud(parula), 'InitialMagnification', 'fit');
                    title(strcat('t = ', num2str(t)));
                    t = t+1;
                end
                h.NextPlot = 'add';
                a = axes; 
                %// Set the title and get the handle to it
                ht = title(strcat('trial ', num2str(trial), ' - layer ', num2str(layer)));
                %// Turn the visibility of the axes off
                a.Visible = 'off';
                %// Turn the visibility of the title on
                ht.Visible = 'on';
            end
        end
        function [X_deci, Y_deci, u_deci, v_deci] = OF(im1, im2, ww)
            w = round(ww/2);

            % Lucas Kanade Here
            % for each point, calculate I_x, I_y, I_t
            Ix_m = conv2(im1,[-1 1; -1 1], 'valid'); % partial on x
            Iy_m = conv2(im1, [-1 -1; 1 1], 'valid'); % partial on y
            It_m = conv2(im1, ones(2), 'valid') + conv2(im2, -ones(2), 'valid'); % partial on t
            u = zeros(size(im1));
            v = zeros(size(im2));

            % within window ww * ww
            for i = w+1:size(Ix_m,1)-w
               for j = w+1:size(Ix_m,2)-w
                  Ix = Ix_m(i-w:i+w, j-w:j+w);
                  Iy = Iy_m(i-w:i+w, j-w:j+w);
                  It = It_m(i-w:i+w, j-w:j+w);

                  Ix = Ix(:);
                  Iy = Iy(:);
                  b = -It(:); % get b here

                  A = [Ix Iy]; % get A here
                  nu = pinv(A)*b; % get velocity here

                  u(i,j)=nu(1);
                  v(i,j)=nu(2);
               end
            end

            % downsize u and v
        %     u_deci = u(1:2:end, 1:2:end);
        %     v_deci = v(1:2:end, 1:2:end);
            u_deci = u;
            v_deci = v;
            % get coordinate for u and v in the original frame
            [m, n] = size(im1);
            [X,Y] = meshgrid(1:n, 1:m);
        %     X_deci = X(1:2:end, 1:2:end);
        %     Y_deci = Y(1:2:end, 1:2:end);
            X_deci = X;
            Y_deci = Y;
        end
        
        function out = histogram2(inputData, nbins)
            out = [];
            minI = min(min(min(inputData)));
            maxI = max(max(max(inputData)));
            for t = 1:size(inputData, 3)
                hist_t_up = []; hist_t_bottom = [];
                for i = 1:(size(inputData, 1)/16)
%                     up = (i*8)-7; center = i*8; bottom = i*16;
                    up = (i*16)-15; center = (i*16)-8; bottom = i*16;
                    tmp = histogram(inputData(up:center,:,t), nbins, 'BinLimits',[minI,maxI]);
                    hist_t_up = cat(1, hist_t_up, tmp.Values);
                    tmp = histogram(inputData(center+1:bottom,:,t), nbins, 'BinLimits',[minI,maxI]);
                    hist_t_bottom = cat(1, hist_t_bottom, tmp.Values);
                end
                out = cat(2, out, hist_t_up, hist_t_bottom);
            end
        end
        
        function out = mean2_sd2_t_col(inputData, issd)
            out = [];
            for t = 1:size(inputData, 3)
                mean_up = []; mean_bottom = []; sd_up = []; sd_bottom = [];
                for i = 1:(size(inputData, 1)/16)
                    up = (i*16)-15; center = (i*16)-8; bottom = i*16;
                    mean_up = cat(1, mean_up, mean(inputData(up:center,:,t)));
                    mean_bottom = cat(1, mean_bottom, mean(inputData(center+1:bottom,:,t)));
                    if issd > 0
                        sd_up = cat(1, sd_up, std(inputData(up:center,:,t)));
                        sd_bottom = cat(1, sd_bottom, std(inputData(center+1:bottom,:,t)));
                    end
                end
                if issd > 0
                    out = cat(2, out, mean_up, mean_bottom, sd_up, sd_bottom);
                else
                    out = cat(2, out, mean_up, mean_bottom);
                end
            end
        end
        
        function out = mean4_sd4_t_col(inputData, issd)
            out = [];
            for t = 1:size(inputData, 3)
                mean_1 = []; mean_2 = []; mean_3 = []; mean_4 = [];
                sd_1 = []; sd_2 = []; sd_3 = []; sd_4 = [];
                for i = 1:(size(inputData, 1)/16)
                    up = (i*16)-15; center = (i*16)-8; bottom = i*16;
                    mean_1 = cat(1, mean_1, mean(inputData(up:center,:,t)));
                    mean_2 = cat(1, mean_2, mean(inputData(center+1:bottom,:,t)));
                    mean_3 = cat(1, mean_3, mean(inputData(up:bottom,1:8,t),2)');
                    mean_4 = cat(1, mean_4, mean(inputData(up:bottom,9:16,t),2)');
                    if issd > 0
                        sd_1 = cat(1, sd_1, std(inputData(up:center,:,t)));
                        sd_2 = cat(1, sd_2, std(inputData(center+1:bottom,:,t)));
                        sd_3 = cat(1, sd_3, std(inputData(up:bottom,1:8,t), [], 2)');
                        sd_4 = cat(1, sd_4, std(inputData(up:bottom,9:16,t), [], 2)');
                    end
                end
                if issd > 0
                    out = cat(2, out, mean_1, mean_2, mean_3, mean_4, sd_1, sd_2, sd_3, sd_4);
                else
                    out = cat(2, out, mean_1, mean_2, mean_3, mean_4);
                end
            end
        end
        
        function varargout = positive_rebound(data, time) % max_mean1
            % pr = positive_rebound(data, time) takes in input the 5D
            % data and the time array and gives the value (in mV) of the
            % highest mean of 16-by-16 matrix in time.
            % [pr, i] = positive_rebound(data, time), gives also the index
            % of the highest mean value.
            
            maximum = []; idx = [];
            
            t_stim_on = find(time==0);
            if isempty(t_stim_on)
               t_stim_on = 1;
               start_t = t_stim_on + 14;
            else
               start_t = t_stim_on + 15;
            end
            
            for c = 1:size(data,1) % c = class
                for i = 1:size(data, 2) % i = sample n.
                    mean1 = featureExtraction.mean_sd_grad_each_t(data(c,i,:, :, start_t:end), 'n_submatrices', 1);
                    [m_tmp, i_tmp] = max(mean1);
                    maximum = cat(1, maximum, m_tmp);
                    idx = cat(1, idx, i_tmp+start_t-1);
                end
            end
            varargout{1} = maximum; varargout{2} = idx;
        end
        
        function varargout = response_peak(data, time) % min_mean1
            minimum = []; idx = [];
            t_stim_on = find(time==0);
            if isempty(t_stim_on)
               t_stim_on = 1;
               start_t = t_stim_on + 3;
            else
               start_t = t_stim_on + 4;
            end
            n_class = size(data,1);
            n_samples = size(data, 2);
            [pr,pr_idx] = featureExtraction.positive_rebound(data, time);
            for c = 1:n_class % c = class
                for i = 1:n_samples
                    mean1 = featureExtraction.mean_sd_grad_each_t(data(c,i,:, :, start_t:pr_idx((c-1)*n_samples + i)), 'n_submatrices', 1);
                    [m_tmp, i_tmp] = min(mean1);
                    minimum = cat(1, minimum, m_tmp);
                    idx = cat(1, idx, i_tmp+start_t-1);
                end
            end
            varargout{1} = minimum; varargout{2} = idx;
        end
        
        function [value, index] = response_onset(data, time, varargin)
            if isempty(varargin)
                t_window = 4;
            elseif length(varargin) == 1
                t_window = varargin{1};
            else
                error('too many inputs')
            end
            t_stim_on = find(time==0);
            if isempty(t_stim_on)
               t_stim_on = 1;
            else
                t_stim_on = t_stim_on +1;
            end
            
            [rp,rp_idx] = featureExtraction.response_peak(data, time);
%             [pr,pr_idx] = featureExtraction.positive_rebound(data, time);
            n_samples = size(data, 2);
            value = []; index = [];
            for c = 1:size(data,1) % c = class
                for i = 1:size(data, 2)
                    mean1 = featureExtraction.mean_sd_grad_each_t(data(c,i,:, :, :), 'n_submatrices', 1);
                    %[response_peak, response_peak_idx] = min(mean1(stim_on+2:end));
                    %response_peak_idx = response_peak_idx+stim_on+1
                    count = 0;
                    l = length(value);
                    [M, M_idx] = max(mean1(t_stim_on:rp_idx((c-1)*n_samples + i)));
                    start_t = t_stim_on + M_idx -1;
                    for t = start_t:rp_idx((c-1)*n_samples + i)-1 %stim_on+3:response_peak_idx
                        %disp(strcat('time point',num2str(t)))
                        if mean1(t+1) < mean1(t)
                            count = count+1;
                            if (count == t_window) % && (mean1(t+2) < mean1(t-t_window+1))
                                value = cat(1, value, mean1(t-t_window+1));
                                index = cat(1, index, t-t_window+1);
                                break
                            end
                        else
                            count = 0;
                        end
                    end
                    if length(value) == l && t_window >= 1 % 1
                        disp(strcat('response onset not found for data #', num2str(i), '. Narrowing the time window to ', num2str(t_window-1)))
                        [tmp_value, tmp_index] = featureExtraction.response_onset(data(c,i,:, :, :), time, t_window-1);
                        value = cat(1, value, tmp_value);
                        index = cat(1, index, tmp_index);
                    elseif length(value) == l && t_window < 1
                        %error(strcat('response onset not found for data #', num2str(i)))
                        value = cat(1, value, mean1(start_t));
                        index = cat(1, index, start_t);
                    end
                end
            end
        end
        
        function out = response_peak_amplitude(data, time, varargin) % response_peak
            if isempty(varargin)
                t_window = 4;
            elseif length(varargin) == 1
                t_window = varargin{1};
            else
                error('too many inputs')
            end
            out = [];
            for c = 1:size(data,1) % c = class
                for i = 1:size(data, 2)
                    response_peak = featureExtraction.response_peak(data(c,i,:, :, :), time);
                    resp_onset = featureExtraction.response_onset(data(c,i,:, :, :), time, t_window);
                    out = cat(1, out, resp_onset-response_peak);
                end
            end
        end

        function out = response_onset_latency(data, time, varargin)
            if isempty(varargin)
                t_window = 4;
            elseif length(varargin) == 1
                t_window = varargin{1};
            else
                error('too many inputs')
            end
            if length(time) ~= size(data, 5)
                error('Time array not consistent with the data')
            end
    % time stim on is 0ms per definition!
%             t_stim_on = find(time==0);
%             if isempty(t_stim_on)
%                error('In time array the stimulus onset is not found')
%             end
            out = [];
            for c = 1:size(data,1) % c = class
                for i = 1:size(data, 2)
                    [resp_onset, time_onset] = featureExtraction.response_onset(data(c,i,:, :, :), time, t_window);
                    out = cat(1, out, time(time_onset));
                end
            end
        end
        
        function out = response_peak_latency(data, time, varargin)
            
            if length(time) ~= size(data, 5)
                error('Time array not consistent with the data')
            end
    % time stim on is 0ms per definition!
%             t_stim_on = find(time==0);
%             if isempty(t_stim_on)
%                error('In time array the stimulus onset is not found')
%             end
            out = [];
            for c = 1:size(data,1) % c = class
                for i = 1:size(data, 2)
                    [resp_peak, resp_peak_idx] = featureExtraction.response_peak(data(c,i,:, :, :), time);
                    out = cat(1, out, time(resp_peak_idx));
                end
            end
        end
        
        function varargout = time_norm_LFP(data, time, varargin)
            if isempty(varargin)
                t_window = 4;
            elseif length(varargin) == 1
                t_window = varargin{1};
            else
                error('too many inputs')
            end
            if length(time) ~= size(data, 5)
                error('Time array not consistent with the data')
            end
            tLFP = []; ctr_resp_onset = []; ctr_time_onset = [];
            for c = 1:size(data,1) % c = class
                for i = 1:size(data, 2)
                    mean1 = featureExtraction.mean_sd_grad_each_t(data(c,i,:, :, :), 'n_submatrices', 1);
                    [resp_onset, time_onset] = featureExtraction.response_onset(data(c,i,:, :, :), time, t_window);
                    t = time_onset + 1;
                    while t < size(data,5) && mean1(t) < mean1(time_onset)
                        t = t+1;
                    end
                    ctr_resp_onset = cat(1, ctr_resp_onset, mean1(t));
                    ctr_time_onset = cat(1, ctr_time_onset, t);
%                     AUC = (time(time_onset+1) - time(time_onset))*trapz(mean1(time_onset:t));
%                     tLFP = cat(1, tLFP, AUC/(time(t)-time(time_onset)));
                    AUC = trapz(mean1(time_onset:t));
                    tLFP = cat(1, tLFP, AUC/(t-time_onset));
                end
            end
            varargout{1} = tLFP;
            varargout{2} = ctr_resp_onset;
            varargout{3} = ctr_time_onset;
        end

    end
end