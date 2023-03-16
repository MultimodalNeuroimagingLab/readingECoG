%% Load data

load( fullfile( dataPath.output, 'preprocessing_0.mat'));

%% ERP - Baseline subtraction

% Removing bad trials (3 in total for subject 1, all during categorization)

% WC 5 - trial 2
data(:, 1, :, 2, 2, 15) = nan;

% WPC 0 - trial 3
data(:, 1, :, 2, 3, 6) = nan;

% WPC 50 - trial 2
data(:, 1, :, 2, 3, 8) = nan;

% new: mean across time -500 to -100 ms, then across 6 trials, all stim, and then 2 tasks
winlen_base = find(epochtime >= -0.5 & epochtime < -0.1);
data_base = mean( mean( mean( mean( data( winlen_base, :, :, :, :, :), 1, 'omitnan'), 5, 'omitnan'), 6, 'omitnan'), 4, 'omitnan');

data_bc = data - data_base;

clear winlen_base data_base

%%  Broadband - Baseline subtraction and normalization

% Removing bad trials (3 in total for subject 1, all during categorization)

% WC 5 - trial 2
bb_data(:, 1, :, 2, 2, 15) = nan;

% WPC 0 - trial 3
bb_data(:, 1, :, 2, 3, 6) = nan;

% WPC 50 - trial 2
bb_data(:, 1, :, 2, 3, 8) = nan;

% new: mean across time -500 to -100 ms, then across 6 trials, all stim, and then 2 tasks
winlen_base = find(epochtime_bb >= -0.5 & epochtime_bb < -0.1);
bb_base = mean( mean( mean( mean( bb_data( winlen_base, :, :, :, :, :), 1, 'omitnan'), 5, 'omitnan'), 6, 'omitnan'), 4, 'omitnan');

% bb x-fold change
bbdata_pc = bsxfun( @rdivide, bb_data, bb_base) - 1;

clear winlen_base bb_base

%%  Alpha - Baseline subtraction and normalization

% Removing bad trials (3 in total for subject 1, all during categorization)

% WC 5 - trial 2
alpha_data(:, 1, :, 2, 2, 15) = nan;

% WPC 0 - trial 3
alpha_data(:, 1, :, 2, 3, 6) = nan;

% WPC 50 - trial 2
alpha_data(:, 1, :, 2, 3, 8) = nan;

% new: mean across time -500 to -100 ms, then across 6 trials, all stim, and then 2 tasks
winlen_base = find(epochtime_bb >= -0.5 & epochtime_bb < -0.1);
alpha_base = mean( mean( mean( mean( alpha_data( winlen_base, :, :, :, :, :), 1, 'omitnan'), 5, 'omitnan'), 6, 'omitnan'), 4, 'omitnan');

% bb x-fold change
alphadata_pc = bsxfun( @rdivide, alpha_data, alpha_base) - 1;

clear winlen_base alpha_base

%% Plot ERP, PSD, BB+scaling, Alpha+scaling

stimgroups  = {[ 6      7      8      9       4  ] [ 10    11     12     13       5  ] [ 14    15     16      4  ] [ 17    18    19      5   ]};    % [20 21 22 10] [2 23 24]};
stimleg     = {["0%", "25%", "50%", "75%", "100%"] ["0%", "25%", "50%", "75%", "100%"] ["3%", "5%",  "8%", "100%"] ["4%", "6%", "10%", "100%"]};    %  'Noise Con'   'Other' };
stimgrnames = {'Word PC'                            'Face PC'                           'Word Con'                  'Face Con'};

% One plot per subject per electrode per stimgroup
for zzz=1%:length(subjects)

    % Ignoring non-recording/non-ECoG/mini channels
    if zzz == 1
        good_channels = [1:26, 29:78];
    elseif zzz == 2
        good_channels = [1:46, 65:110];
    else
        assert(1==0);
    end
   
    for ccc = good_channels
        % Rest PSD: Avg across trials (6), then stim (18), then tasks (2)
        psd_rest = squeeze( mean( mean( mean( psd_off( zzz, ccc, :, :, :, :), 4, 'omitnan'), 5, 'omitnan'), 3, 'omitnan'));
        %psd_rest = movmean( psd_rest, 10, 1);
            
        for stimg = 1:length(stimgroups)
            
            figure('Position', [0 0 250*length(stimgroups{stimg}) 1500], 'Visible', 'off')
            
            for stimc = 1:length(stimgroups{stimg})

                %%%%%% ERP

                subplot(6, length(stimgroups{stimg}), stimc);
                hold on;

                winlen_tc = find(epochtime >= -0.5 & epochtime <= 2.5);   % -0.5 s to 2.5 s

                % Fixation - normalizing and smoothing 10x to 200Hz/5ms resolution
                data_smooth = movmean( squeeze( data_bc( winlen_tc, zzz, ccc, 1, :, stimgroups{stimg}(stimc))), 1, 1, "omitnan");

                nTrials = size( data_smooth, 2) - sum( all( isnan(data_smooth), 1));
        
                yMean = mean( data_smooth, 2, 'omitnan');
                ySEM  = std( data_smooth, 0, 2, 'omitnan') / (sqrt(nTrials));
                ts = tinv([0.025, 0.975], nTrials - 1);         % n-1 degrees of freedom -HH
                conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
                plot(yMean, 'b-', 'LineWidth', 2);
                xline(find(epochtime == 0), '--k')
                fill([1:length(yMean) length(yMean):-1:1],conf, [0 0 .5],'EdgeColor','none','FaceAlpha',.3,'HandleVisibility','off');

                % Categorization - normalizing and smoothing 10x to 200Hz/5ms resolution
                data_smooth = movmean( squeeze( data_bc( winlen_tc, zzz, ccc, 2, :, stimgroups{stimg}(stimc))), 1, 1, "omitnan");

                nTrials = size( data_smooth, 2) - sum( all( isnan(data_smooth), 1));
                
                yMean = mean( data_smooth, 2, 'omitnan');
                ySEM  = std( data_smooth, 0, 2, 'omitnan') / (sqrt(nTrials));
                ts = tinv([0.025, 0.975], nTrials - 1);         % n-1 degrees of freedom -HH
                conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
                plot(yMean, 'r-', 'LineWidth', 2);
                xline(find(epochtime == 0), '--k')
                fill([1:length(yMean) length(yMean):-1:1],conf, [.5 0 0],'EdgeColor','none','FaceAlpha',.3,'HandleVisibility','off');
                ylabel('Normalized Amplitude');
                xlabel('Time (s)');
                set(gca,'XTick',1:1000:length(winlen_tc)); set(gca,'XTickLabel',epochtime(winlen_tc(1)):0.5:epochtime(winlen_tc(end)));
                title(stimgrnames{stimg}+" "+stimleg{stimg}(stimc)); 

                hold off

                %%%%%% PSD

                subplot(6, length(stimgroups{stimg}), length(stimgroups{stimg})+stimc);
                hold on;

                % Rest
                plot(f_psd, normalize( 10*log10(psd_rest)), 'k-', 'LineWidth', 2);

                % Fixation
                temp_psd = squeeze( psd_on( zzz, ccc, 1, :, stimgroups{stimg}(stimc), :))';
                %temp_psd = movmean( temp_psd, 10, 1);

                yMean = mean( temp_psd, 2, 'omitnan');
                plot(f_psd, normalize( 10*log10(yMean)), 'b-', 'LineWidth', 2);
                
                % Categorization
                temp_psd = squeeze( psd_on( zzz, ccc, 2, :, stimgroups{stimg}(stimc), :))';
                %temp_psd = movmean( temp_psd, 10, 1);

                yMean = mean( temp_psd, 2, 'omitnan');
                plot(f_psd, normalize( 10*log10(yMean)), 'r-', 'LineWidth', 2);

                ylabel('Normalized PSD');
                xlabel('Frequency (Hz)');

                hold off

            end

            %%%%%% BB Timeseries

            mx = -Inf; mn = Inf;

            for stimc = 1:length(stimgroups{stimg})

                subplot(6, length(stimgroups{stimg}), 2*length(stimgroups{stimg})+stimc);
                hold on;

                winlen_tc = find(epochtime_bb >= -0.5 & epochtime_bb <= 2.5);   % -0.5 s to 2.5 s
                
                % Fixation - smoothing 10x to 20Hz/50ms resolution
                data_smooth = movmean( squeeze( bbdata_pc( winlen_tc, zzz, ccc, 1, :, stimgroups{stimg}(stimc))), 10, 1, "omitnan");

                nTrials = size( data_smooth, 2) - sum( all( isnan(data_smooth), 1));
        
                yMean = mean( data_smooth, 2, 'omitnan');
                ySEM  = std( data_smooth, 0, 2, 'omitnan') / (sqrt(nTrials));
                ts = tinv([0.025, 0.975], nTrials - 1);         % n-1 degrees of freedom -HH
                conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
                plot(yMean, 'b-', 'LineWidth', 2);
                fill([1:length(yMean) length(yMean):-1:1],conf, [0 0 .5],'EdgeColor','none','FaceAlpha',.3,'HandleVisibility','off');

                % Categorization - smoothing 10x to 20Hz/50ms resolution
                data_smooth = movmean( squeeze( bbdata_pc( winlen_tc, zzz, ccc, 2, :, stimgroups{stimg}(stimc))), 10, 1, "omitnan");

                nTrials = size( data_smooth, 2) - sum( all( isnan(data_smooth), 1));
                
                yMean = mean( data_smooth, 2, 'omitnan');
                ySEM  = std( data_smooth, 0, 2, 'omitnan') / (sqrt(nTrials));
                ts = tinv([0.025, 0.975], nTrials - 1);         % n-1 degrees of freedom -HH
                conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
                plot(yMean, 'r-', 'LineWidth', 2);
                xline(find(epochtime_bb == 0), '--k')
                fill([1:length(yMean) length(yMean):-1:1],conf, [.5 0 0],'EdgeColor','none','FaceAlpha',.3,'HandleVisibility','off');
                ylabel('BB (x-fold increase)');
                xlabel('Time (s)');
                set(gca,'XTick',1:100:length(winlen_tc)); set(gca,'XTickLabel',epochtime_bb(winlen_tc(1)):0.5:epochtime_bb(winlen_tc(end)));

                ax = axis;
                mx = max(mx,max(ax(3:4)));
                mn = min(mn,min(ax(3:4)));

                hold off

            end

            for stimc = 1:length(stimgroups{stimg})

                subplot(6, length(stimgroups{stimg}), 2*length(stimgroups{stimg})+stimc);
                ax = axis;
                axis([ax(1:2) mn mx]);

            end

            %%%%%% BB Scaling Factor

            mx = -Inf; mn = Inf;

            for stimc = 1:length(stimgroups{stimg})

                subplot(6, length(stimgroups{stimg}), 3*length(stimgroups{stimg})+stimc);
                hold on;

                winlen_sf = find(epochtime_bb >= 0 & epochtime_bb <= 2);        %    0 s to 1 s

                scaling_factor = movmean( ...
                                    zerodiv( ...
                                        squeeze( bbdata_pc( winlen_sf, zzz, ccc, 2, :, stimgroups{stimg}(stimc))) +1, ...
                                        squeeze( bbdata_pc( winlen_sf, zzz, ccc, 1, :, stimgroups{stimg}(stimc))) +1, ...
                                        1), ...
                                    10, 1, "omitnan");

                nTrials = size(scaling_factor, 2) - sum( all( isnan(scaling_factor), 1));
                yMean = mean(scaling_factor, 2, "omitnan");
                ySEM  = std(scaling_factor, 0, 2, "omitnan") / sqrt(nTrials);
                ts = tinv([0.025, 0.975], nTrials - 1); % n-1 degrees of freedom -HH
                conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
                plot(yMean,'k-', 'LineWidth', 2);
                fill([1:length(winlen_sf) length(winlen_sf):-1:1],conf, [.5 0.5 0.5],'EdgeColor','none','FaceAlpha',.5,'HandleVisibility','off');
                yline(1, '--k')
                ylabel('BB Scaling (x-fold)');
                xlabel('Time (s)');
                set(gca,'XTick',1:100:length(winlen_sf)); set(gca,'XTickLabel',epochtime_bb(winlen_sf(1)):0.5:epochtime_bb(winlen_sf(end)));

                ax = axis;
                mx = max(mx,max(ax(3:4)));
                mn = min(mn,min(ax(3:4)));

                hold off

            end

            for stimc = 1:length(stimgroups{stimg})

                subplot(6, length(stimgroups{stimg}), 3*length(stimgroups{stimg})+stimc);
                ax = axis;
                axis([ax(1:2) mn mx]);

            end

            %%%%%% Alpha Timeseries

            mx = -Inf; mn = Inf;

            for stimc = 1:length(stimgroups{stimg})

                subplot(6, length(stimgroups{stimg}), 4*length(stimgroups{stimg})+stimc);
                hold on;

                winlen_tc = find(epochtime_bb >= -0.5 & epochtime_bb <= 2.5);   % -0.5 s to 2.5 s
                
                % Fixation - smoothing 10x to 20Hz/50ms resolution
                data_smooth = movmean( squeeze( alphadata_pc( winlen_tc, zzz, ccc, 1, :, stimgroups{stimg}(stimc))), 10, 1, "omitnan");

                nTrials = size( data_smooth, 2) - sum( all( isnan(data_smooth), 1));
        
                yMean = mean( data_smooth, 2, 'omitnan');
                ySEM  = std( data_smooth, 0, 2, 'omitnan') / (sqrt(nTrials));
                ts = tinv([0.025, 0.975], nTrials - 1);         % n-1 degrees of freedom -HH
                conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
                plot(yMean, 'b-', 'LineWidth', 2);
                fill([1:length(yMean) length(yMean):-1:1],conf, [0 0 .5],'EdgeColor','none','FaceAlpha',.3,'HandleVisibility','off');

                % Categorization - smoothing 10x to 20Hz/50ms resolution
                data_smooth = movmean( squeeze( alphadata_pc( winlen_tc, zzz, ccc, 2, :, stimgroups{stimg}(stimc))), 10, 1, "omitnan");

                nTrials = size( data_smooth, 2) - sum( all( isnan(data_smooth), 1));
                
                yMean = mean( data_smooth, 2, 'omitnan');
                ySEM  = std( data_smooth, 0, 2, 'omitnan') / (sqrt(nTrials));
                ts = tinv([0.025, 0.975], nTrials - 1);         % n-1 degrees of freedom -HH
                conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
                plot(yMean, 'r-', 'LineWidth', 2);
                xline(find(epochtime_bb == 0), '--k')
                yline(0, '--k')
                fill([1:length(yMean) length(yMean):-1:1],conf, [.5 0 0],'EdgeColor','none','FaceAlpha',.3,'HandleVisibility','off');
                ylabel('Alpha (x-fold increase)');
                xlabel('Time (s)');
                set(gca,'XTick',1:100:length(winlen_tc)); set(gca,'XTickLabel',epochtime_bb(winlen_tc(1)):0.5:epochtime_bb(winlen_tc(end)));
                title(stimgrnames{stimg}+" "+stimleg{stimg}(stimc)); 

                ax = axis;
                mx = max(mx,max(ax(3:4)));
                mn = min(mn,min(ax(3:4)));

                hold off

            end

            for stimc = 1:length(stimgroups{stimg})

                subplot(6, length(stimgroups{stimg}), 4*length(stimgroups{stimg})+stimc);
                ax = axis;
                axis([ax(1:2) mn mx]);

            end
                
            %%%%%% Alpha Scaling Factor
            
            mx = -Inf; mn = Inf;

            for stimc = 1:length(stimgroups{stimg})

                subplot(6, length(stimgroups{stimg}), 5*length(stimgroups{stimg})+stimc);
                hold on;

                winlen_sf = find(epochtime_bb >= 0 & epochtime_bb <= 2);        %    0 s to 1 s

                scaling_factor = movmean( ...
                                    zerodiv( ...
                                        squeeze( alphadata_pc( winlen_sf, zzz, ccc, 2, :, stimgroups{stimg}(stimc))) +1, ...
                                        squeeze( alphadata_pc( winlen_sf, zzz, ccc, 1, :, stimgroups{stimg}(stimc))) +1, ...
                                        1), ...
                                    10, 1, "omitnan");

                nTrials = size(scaling_factor, 2) - sum( all( isnan(scaling_factor), 1));
                yMean = mean(scaling_factor, 2, "omitnan");
                ySEM  = std(scaling_factor, 0, 2, "omitnan") / sqrt(nTrials);
                ts = tinv([0.025, 0.975], nTrials - 1); % n-1 degrees of freedom -HH
                conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
                plot(yMean,'k-', 'LineWidth', 2);
                fill([1:length(winlen_sf) length(winlen_sf):-1:1],conf, [.5 0.5 0.5],'EdgeColor',[.5 0.5 0.5],'FaceAlpha',.5,'HandleVisibility','off');
                yline(1, '--k')
                ylabel('Alpha Scaling (x-fold)');
                xlabel('Time (s)');
                set(gca,'XTick',1:100:length(winlen_sf)); set(gca,'XTickLabel',epochtime_bb(winlen_sf(1)):0.5:epochtime_bb(winlen_sf(end)));

                ax = axis;
                mx = max(mx,max(ax(3:4)));
                mn = min(mn,min(ax(3:4)));

                hold off

            end

            for stimc = 1:length(stimgroups{stimg})

                subplot(6, length(stimgroups{stimg}), 5*length(stimgroups{stimg})+stimc);
                ax = axis;
                axis([ax(1:2) -10 10]);

            end

            figurewrite(sprintf('sub%d_ch%03d_%s', zzz, ccc, stimgrnames{stimg}),[],[],fullfile(dataPath.output, 'timeseries', sprintf('sub%d_ch%03d', zzz, ccc)));
            close all
        end
    end
end