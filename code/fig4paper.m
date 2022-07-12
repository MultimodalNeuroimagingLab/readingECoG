%% BB x-fold change for EVC (100% con/phase)

% Do we want fixation and categorization together or separate?

stimgroups  = {[5], [4]};
stimleg  = {"Face", "Word"};

elec_consider{1} = [77];
elec_consider{2} = [110];
winlen = 1:601; % -0.5 to 2.5 s

for zzz = 1:2
    for ccc = elec_consider{zzz}
        figure; hold on;
        for ppp = 1:2
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

%% Plotting Spectra

for zzz = 1:2
    for ccc = elec_consider{zzz}
        for ppp = 1:2
            for stimg = 1:length(stimgroups)
                for stimc = 1:length(stimgroups{stimg})
                    temp = mean(data(:, zzz, ccc, ppp, :, stimgroups{stimg}(stimc)), 5);
                    temp(isnan(temp)) = 0;
                    [S,f_spectra] = getWaveletSpectrogram(temp, fsorig, [1, fupper]);
                    spectra = log10(S) - log10(squeeze(spect_base(zzz,ccc,:)));
                    %plotting spectra
                    load('dg_colormap.mat', 'cm') %KM color scheme;
                    figure, uimagesc(epochtime(1:9001),f_spectra,spectra,[-1 1]); axis xy; colormap(cm); colorbar
                    xlabel('Time (s)');
                    ylabel('Frequency (Hz)');
                    title( stimgrnames{stimg} + " " + stimleg{stimg}(stimc))
                    %figurewrite(stimgrnames{stimg} + " " + stimleg{stimg}(stimc), -1, [], sprintf('spectrogram/Subj%d_ch%03d_%s_%s', zzz, ccc, channellabels{zzz}{ccc}, tasklabels{ppp}));
                end
            end
        end
    end
end

%% Fixation v/s Categorization over con and phase coherence averaged over all visually responsive electrodes and face and word category (general overview)

%gcc = {horzcat(9, 15, 41, 44, 74, 76), horzcat(21, 44, 87, 110) };
%non-selective electrodes

stimgroups  = {[6 7 8 9 4] [10 11 12 13 5] [14 15 16 4] [17 18 19 5]};
stimleg  = {["0%", "25%", "50%", "75%", "100%"] ["4%", "6%", "10%", "100%"]};
stimgrnames = {'Phase Coherence' 'Contrast'};

winlen = 1:601; % -0.5 to 2.5 s

fix_avgresp = zeros(length(epochtime_bb), length(gcc{1})+length(gcc{2}), 2*numtrials, 2, 5); % fixation - winlen x number of visually responsive electrodes x 2*number of trials x #of stim in each category
cat_avgresp = zeros(length(epochtime_bb), length(gcc{1})+length(gcc{2}), 2*numtrials, 2, 5); % categorization

% Combining over subjects and electrodes
counter = 1;
for zzz = 1:2
    for ccc = gcc{zzz}
        for stimg = 1:2 % we need only 1st and 3rd group
            for stimc = 1 : length(stimgroups{2*stimg-1})
                fix_avgresp(:, counter,  1:6, stimg, stimc) = squeeze( bbdata_pc(:, zzz, ccc, 1, :, stimgroups{2*stimg-1}(stimc))); % fixation 
                cat_avgresp(:, counter, 7:12, stimg, stimc) = squeeze( bbdata_pc(:, zzz, ccc, 2, :, stimgroups{2*stimg-1}(stimc))); % categorization
            end
        end
        counter = counter+1;
    end
%  temp3(:,:,zzz) = squeeze( mean( mean( bbdata_pc(winlen, zzz, gcc{zzz}, 1, :, :), 5) ,3));
end

fix_avgresp = reshape(fix_avgresp, [length(epochtime_bb), (length(gcc{1})+length(gcc{2}))*2*numtrials, 2, 5]);
cat_avgresp = reshape(cat_avgresp, [length(epochtime_bb), (length(gcc{1})+length(gcc{2}))*2*numtrials, 2, 5]);

% Plot


for stimg = 1:2
    
    figureprep([200 200 3800 2200]);
    mx = -Inf; mn = Inf;
    
    for stimc = 1 : length(stimgroups{2*stimg-1})
        
        subplot(2,length(stimgroups{2*stimg-1}),stimc); hold on;
        
        nTrials = size(fix_avgresp(:,:,stimg,stimc),2);
        
        % Fixation
        yMean = smooth(mean(fix_avgresp(:,:,stimg,stimc),2),40);
        ySEM  = smooth(std(fix_avgresp(:,:,stimg,stimc), 0, 2)/(6*sqrt(nTrials)),40);
        ts = tinv([0.025, 0.975], nTrials - 1); % n-1 degrees of freedom -HH
        conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
        plot(smooth(yMean,40),'b-', 'LineWidth', 2);
        fill([1:length(epochtime_bb) length(epochtime_bb):-1:1],conf, [0 0 .5],'EdgeColor',[0 0 .5],'FaceAlpha',.5,'HandleVisibility','off');
        
        % Categorization
        yMean = smooth(mean(cat_avgresp(:,:,stimg,stimc),2),40);
        ySEM  = smooth(std(cat_avgresp(:,:,stimg,stimc), 0, 2)/(6*sqrt(nTrials)),40);
        ts = tinv([0.025, 0.975], nTrials - 1); % n-1 degrees of freedom -HH
        conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
        plot(smooth(yMean,40),'r-', 'LineWidth', 2);
        fill([1:length(epochtime_bb) length(epochtime_bb):-1:1],conf, [.5 0 0],'EdgeColor',[.5 0 0],'FaceAlpha',.5,'HandleVisibility','off');
        
        ax = axis;
        mx = max(mx,max(ax(3:4)));
        mn = min(mn,min(ax(3:4)));
        xlim([0 800]);
        ylabel('BB response (x-fold)');
        xlabel('t (s)');
        set(gca,'XTick',0:100:winlen(end)); set(gca,'XTickLabel',-0.5:0.5:epochtime_bb(winlen(end)));
        title(stimgrnames{stimg}+" "+stimleg{stimg}(stimc));
        legend('F','C');
        
    end
    
    for stimc = 1:length(stimgroups{2*stimg-1})
        
        subplot(2,length(stimgroups{2*stimg-1}),stimc); hold on;
        ax = axis;
        axis([ax(1:2) mn mx]);
        %exportgraphics(gca,[sprintf('FvsC/%s%s/',stimgrnames{stimg},stimleg{stimg}(stimc)) 'STC.eps']);
    
    end
        
    mx = -Inf; mn = Inf;
    
    for stimc = 1 : length(stimgroups{2*stimg-1})
        
        subplot(2,length(stimgroups{2*stimg-1}), length(stimgroups{2*stimg-1})+stimc); hold on;
        
        scaling_factor = squeeze(cat_avgresp(:,:,stimg,stimc)) - squeeze(fix_avgresp(:,:,stimg,stimc));
        
        nTrials = size(scaling_factor,2);
        yMean = smooth(mean(scaling_factor,2),40);
        ySEM  = smooth(std(scaling_factor,0,2)/(6*sqrt(nTrials)),40);
        ts = tinv([0.025, 0.975], nTrials - 1); % n-1 degrees of freedom -HH
        conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
        plot(yMean,'k-', 'LineWidth', 2);
        fill([1:length(yMean) length(yMean):-1:1],conf, [.5 0.5 0.5],'EdgeColor',[.5 0.5 0.5],'FaceAlpha',.5,'HandleVisibility','off');
        %[h,p,ci,stats] = ttest2(mean(stimresp_c(:,:,2,stimg,stimc),2),mean(avgresp(:,:,2,stimg,stimc),2));
        for t = 1:length(yMean)     % CoCoSys Lab
            [~, pval(t)] = ttest(squeeze(cat_avgresp(t,:,stimg,stimc)), squeeze(fix_avgresp(t,:,stimg,stimc)));
        end
        % convert to logical
        signific = nan(1, length(yMean)); 
        signific(pval < 0.05) = 1;
        plot(signific * -0.03, '.k');
        % indicate what we're showing
        %text(10.2, -0.03, 'p < 0.001');
        ax = axis;
        axis([100 300 ax(3:4)]);
        mx = max(mx,max(ax(3:4)));
        mn = min(mn,min(ax(3:4)));
        %xlim([100 300]);
        ylabel('BB response (x-fold)');
        xlabel('t (s)');
        set(gca,'XTick',0:100:400); set(gca,'XTickLabel',-0.5:0.5:1.5);
        %title(stimgrnames{stimg}+" "+stimleg{stimg}(stimc));
        title('Scaling Factor');
        
    end
    
    for stimc = 1:length(stimgroups{2*stimg-1})
        
        subplot(2,length(stimgroups{2*stimg-1}), length(stimgroups{2*stimg-1})+stimc); hold on;
        ax = axis;
        axis([ax(1:2) mn mx]);
        %exportgraphics(gca,[sprintf('FvsC/%s%s/',stimgrnames{stimg},stimleg{stimg}(stimc)) 'SF.eps']);

    end
    figurewrite(sprintf('%s',stimgrnames{stimg}),[],[],'FvsC');
end

