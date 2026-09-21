%% Preprocessing_1 
% 
% Load curated data from derivatives, remove bad trials, baseline subtraction and normalization
%
% ZQ 2023 Mayo Clinic
%
%% Load curated data
%
% for details on whats inside data_ecogTaskDemands.mat, refer the readme document

close all; clear; clc; 

% set your data path
dataPath = setLocalDataPath(1); 
load( fullfile( dataPath.input, 'derivatives', 'data_ecogTaskDemands.mat'));

%% ERP - Baseline subtraction

% Removing bad trials (4 in total for subject 1,  3 of 4 during categorization)

% WC 5 - trial 2 - cat
data(:, 1, :, 2, 2, 15) = nan;

% WPC 0 - trial 3 - cat
data(:, 1, :, 2, 3, 6) = nan;

% WPC 50 - trial 3 - cat
data(:, 1, :, 2, 3, 8) = nan;

% WC/WPC 100 - trial 5 - fix
data(:, 1, :, 1, 5, 4) = nan;

% baseline window: median across time -500 to -100 ms, then across 6 trials, 22 stim, and then 2 tasks
winlen_base = find(epochtime >= -0.5 & epochtime < -0.1);
data_base = median( median( median( median( data( winlen_base, :, :, :, :, :), 1, 'omitnan'), 5, 'omitnan'), 6, 'omitnan'), 4, 'omitnan');

% now baseline subtract
data_bc = data - data_base;

clear winlen_base data_base

%%  HFB - Baseline subtraction and normalization

% Removing bad trials (4 in total for subject 1, 3 of 4 during categorization)

% WC 5 - trial 2
hfb_data(:, 1, :, 2, 2, 15) = nan;

% WPC 0 - trial 3
hfb_data(:, 1, :, 2, 3, 6) = nan;

% WPC 50 - trial 3
hfb_data(:, 1, :, 2, 3, 8) = nan;

% WC/WPC 100 - trial 5 - fix
hfb_data(:, 1, :, 1, 5, 4) = nan;

% % across baseline (-500 to -100 ms), merge task x trials x stim
winlen_base = find(epochtime_bb >= -0.5 & epochtime_bb < -0.1);
temp_base = reshape( hfb_data( winlen_base, :, :, :, :, :), length(winlen_base), size(hfb_data,2), size(hfb_data,3), size(hfb_data,4) * size(hfb_data,5) * size(hfb_data,6));

% then geomean across the merged trials and time per subject per channel
hfb_base = geomean( geomean( temp_base, 4, 'omitnan'), 1, 'omitnan');

% hfb x-fold change
% eg. HFB power of 8 wrt a baseline of 3 is now 8/3-1 = 2-fold change wrt 0
hfbdata_pc = bsxfun( @rdivide, hfb_data, hfb_base) - 1;

clear winlen_base hfb_base

%%  LFB (alpha+beta) - Baseline subtraction and normalization
% similar to HFB (above)

% Removing bad trials (4 in total for subject 1, 3 of 4 during categorization)

% WC 5 - trial 2
lfb_data(:, 1, :, 2, 2, 15) = nan;

% WPC 0 - trial 3
lfb_data(:, 1, :, 2, 3, 6) = nan;

% WPC 50 - trial 3
lfb_data(:, 1, :, 2, 3, 8) = nan;

% WC/WPC 100 - trial 5 - fix
lfb_data(:, 1, :, 1, 5, 4) = nan;

winlen_base = find(epochtime_bb >= -0.5 & epochtime_bb < -0.1);
temp_base = reshape( lfb_data( winlen_base, :, :, :, :, :), length(winlen_base), size(lfb_data,2), size(lfb_data,3), size(lfb_data,4) * size(lfb_data,5) * size(lfb_data,6));
lfb_base = geomean( geomean( temp_base, 4, 'omitnan'), 1, 'omitnan');

% lfb x-fold change
lfbdata_pc = bsxfun( @rdivide, lfb_data, lfb_base) - 1;

clear winlen_base lfb_base

%% Summary plots - for quick look of whats in the data
%
% feed your curosity!!

% sub-ID 1 or 2
temp_sub = 1;

% task # 1=F, 2=C
temp_task = 2;

if temp_task == 1, temp_tskName = 'Fixation'; 
elseif temp_task == 2, temp_tskName = 'Categorization'; end

% reshape HFB reasponses to time x channel x [stimulus x trials]
temp_sig = squeeze(reshape(hfbdata_pc(:, temp_sub, :, temp_task, :, :), size(hfbdata_pc,1), size(hfbdata_pc,3), []));

% geomean across [stimulus x trials]

temp_sig = squeeze(geomean(temp_sig+1, 3, 'omitnan')-1);

figure
imagesc( [], [], temp_sig',[-1 4])
colormap(cmapsign4)

% set axis for interperabilty
temp_epoch = 1:50:length(epochtime_bb);
set(gca,'XTick',temp_epoch); set(gca,'XTickLabel',epochtime_bb(temp_epoch));
xlabel('Time (s)', 'FontSize',18,'FontWeight','bold')
xline(find(epochtime_bb == 0), 'w', 'LineWidth',2);
xline(find(epochtime_bb == 2), 'w', 'LineWidth',2);
xline(find(epochtime_bb == 4), 'w', 'LineWidth',2);

set(gca,'yTick',1:128); set(gca,'yTickLabel',channellabels{temp_sub})
ylabel('Channel Names', 'FontSize',18,'FontWeight','bold')

colorbar
ylabel(colorbar,'BB Power','FontSize',12,'Rotation',270)
title(['Sub - ', num2str(temp_sub), ' ' temp_tskName, ' BB Power (x-fold change)'])

clear temp*


%% Individual channel traces - quick look

% sub-ID 1 or 2
temp_sub = 1;

% task # 1=F, 2=C
temp_task = 2;

% list of channels you want to see for the above subject & task
vis_el = 77;

% epoch times for ERP
winlen_tt = find(epochtime >= -0.5 & epochtime <= 2.5);   % -0.5 s to 2.5 s

% epoch times for HFB
winlen_tc = find(epochtime_bb >= -0.5 & epochtime_bb <= 2.5);   % -0.5 s to 2.5 s

for ii = vis_el

    % ERP

    figure

    subplot(2,1,1)
    hold on
    
    data_smooth = reshape( data_bc( winlen_tt, temp_sub, ii, temp_task, :, :), length(winlen_tt), []);
    
    nTrials = size( data_smooth, 2) - sum( all( isnan(data_smooth), 1));
    
    % Mean, SEM and 95% CI
    yMean = mean( data_smooth, 2, "omitnan");
    ySEM  = std( data_smooth, 0, 2, "omitnan") / sqrt(nTrials);
    ts = tinv([0.025, 0.975], nTrials - 1);
    conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];

    plot(epochtime(winlen_tt), yMean);
    fill([epochtime(winlen_tt) epochtime(winlen_tt(end:-1:1))],conf, [0 0 0.5],'EdgeColor','none','FaceAlpha',.2,'HandleVisibility','off');

    xlabel('Time (s)')
    ylabel('ERP (uV)')
    title(channellabels{temp_sub}{ii})
    yline(0, '--k')
    xline(0, '--k')
    xline(2, '--k')


    % BB

    subplot(2,1,2)
    hold on
    
    data_smooth = reshape( hfbdata_pc( winlen_tc, temp_sub, ii, temp_task, :, :), length(winlen_tc), []);

    % log-transform as HFB responses typically log-normal distribution
    data_smooth = log(data_smooth+1);
    
    nTrials = size( data_smooth, 2) - sum( all( isnan(data_smooth), 1));
            
    % Mean, SEM and 95% CI
    yMean = mean( data_smooth, 2, "omitnan");
    ySEM  = std( data_smooth, 0, 2, "omitnan") / sqrt(nTrials);
    ts = tinv([0.025, 0.975], nTrials - 1);
    conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];

    % exp to plot in linear scale
    plot(epochtime_bb(winlen_tc), exp(yMean)-1);
    fill([epochtime_bb(winlen_tc) epochtime_bb(winlen_tc(end:-1:1))],exp(conf)-1, [0 0 0.5],'EdgeColor','none','FaceAlpha',.2,'HandleVisibility','off');    

    xlabel('Time (s)')
    ylabel('BB (x-fold increase)')
    yline(0, '--k')
    xline(0, '--k')
    xline(2, '--k')


end

clear temp* ii vis_el data_smooth yMean ySEM ts conf nTrials

