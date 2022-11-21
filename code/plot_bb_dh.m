%%  Baseline subtraction and normalization

% load broadband power
% save('Feb21_exttime_psd.mat','-v7.3');
% load('Feb21_exttime_psd.mat');

% new: mean across time -500 to 0 ms, then across 6 trials, all stim, and then 2 tasks
time_baseline = find(epochtime_bb<0); % -500:0 ms
% average baseline across time, trials, stimuli and tasks, gives 1 value
% per subject per electrode
bb_base = nanmean( nanmean( nanmean( nanmean( bb_data(time_baseline, :, :, :, :, :), 1), 5), 6), 4);
% divide by baseline
bbdata_pc = bsxfun(@rdivide, bb_data, bb_base) - 1;
% --> now bbdata_pc is in percent signal change of the bb power wrt
% baseline

%% timeseries for one electrode 

% stimgroups ={[6     7     8     9]    [10    11    12     13] [14   15  16  4]    [17  18   19     5]};   % [20 21 22 10] [2 23 24]};
% stimleg  = {["0", "25", "50", "75"]  ["0", "25", "50", "75"] ["3" "5" "8" "100"] ["4" "6" "10" "100"]};
% stimgrnames = {'Word Phase'           'Face Phase'           'Word Contrast'     'Face Contrast'};%   'Noise Con'   'Other'};

stimgroups  = {[5], [4]};
stimleg  = {"Face", "Word"};

elec_consider{1} = [77];
elec_consider{2} = [110];

for zzz = 1%:2 % subject
    for ccc = elec_consider{zzz} % electrode
        figure; hold on;
        for ppp = 1:2 % task
            sig_consider = movmean(squeeze(bbdata_pc(:, zzz, ccc, ppp, :, stimgroups{zzz})), 40, 1);
            nTrials = size(sig_consider, 2);
            yMean = mean(sig_consider, 2);
            ySEM  = std(sig_consider, 0, 2) / sqrt(nTrials);
            ts = tinv([0.025, 0.975], nTrials - 1); % n-1 degrees of freedom -HH
            conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
            if ppp == 1
                plot(yMean,'b-', 'LineWidth', 2);
                fill([1:length(epochtime_bb) length(epochtime_bb):-1:1],conf, [0 0 .5],'EdgeColor',[0 0 .5],'FaceAlpha',.5,'HandleVisibility','off');
            elseif ppp == 2
                plot(yMean,'r-', 'LineWidth', 2);
                fill([1:length(epochtime_bb) length(epochtime_bb):-1:1],conf, [.5 0 0],'EdgeColor',[.5 0 0],'FaceAlpha',.5,'HandleVisibility','off');
            end
            
        end
    end
end


%% average responses for all conditions 

% stimgroups ={[6     7     8     9]    [10    11    12     13] [14   15  16  4]    [17  18   19     5]};   % [20 21 22 10] [2 23 24]};
% stimleg  = {["0", "25", "50", "75"]  ["0", "25", "50", "75"] ["3" "5" "8" "100"] ["4" "6" "10" "100"]};
% stimgrnames = {'Word Phase'           'Face Phase'           'Word Contrast'     'Face Contrast'};%   'Noise Con'   'Other'};

plotting_order =[6     7     8     9 10    11    12     13 14   15  16  4 17  18   19     5];   % [20 21 22 10] [2 23 24]};

elec_consider{2} = [86];
winlen = find(epochtime_bb>0 & epochtime_bb<.2);

for zzz = 2%:2 % subject
    for ccc = elec_consider{zzz} % electrode
        figure; hold on;
        % find plotting min/max here
        sig_minmax = squeeze(mean(bbdata_pc(winlen, zzz, ccc, :, :, plotting_order),1));
        yMean = mean(sig_consider, 1);
        axes_minmax = [min(yMean(:)) 1.2*max(yMean(:))];
        for ppp = 1:2 % task
            sig_consider = squeeze(mean(bbdata_pc(winlen, zzz, ccc, ppp, :, plotting_order),1));
            nTrials = size(sig_consider, 1);
            yMean = mean(sig_consider, 1);
            ySEM  = std(sig_consider, 0, 1) / sqrt(nTrials);
            ts = tinv([0.025, 0.975], nTrials - 1); % n-1 degrees of freedom -HH
            conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
            if ppp == 1
                subplot(3,1,1),hold on
                bar(yMean,'w');
                errorbar(1:16,yMean,2*ySEM','k.');
                plot([4.5 4.5],[0 5],'b','LineWidth',2) % stimulus separations
                ylim([0 axes_minmax(2)])
            elseif ppp == 2
                subplot(3,1,2),hold on
                bar(yMean,'w');
                errorbar(1:16,yMean,2*ySEM','k.');
                plot([4.5 4.5],[0 5],'b','LineWidth',2) % stimulus separations
                ylim([0 axes_minmax(2)])
            end
        end
        sig_consider1 = squeeze(mean(bb_data(winlen, zzz, ccc, 1, :, plotting_order),1));
        sig_consider2 = squeeze(mean(bb_data(winlen, zzz, ccc, 2, :, plotting_order),1));
        scaling_factor = mean(sig_consider2, 1)./mean(sig_consider1, 1) ;
        subplot(3,1,3),hold on
        bar(scaling_factor,'w');   
        plot([0 17],[1 1],'k')
        plot([4.5 4.5; 8.5 8.5]',[0 2;0 2 ]','b','LineWidth',2) % stimulus separations
    end
end
