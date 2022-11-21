%% One example for each subject along with spectra

% Do we want both fixation and categorization together or just one ?

% Example electrode selected to show raw-voltage traces 
elec_consider{1} = [46];     % EVC: 75 , VTC: 46
elec_consider{2} = [110];    % LVC: 108, VTC: 87


%% Computing Spectra baseline

winlen_base = find(epochtime > -0.5 & epochtime < -0.1);

stimgroups  = {[6 7 9]   [10 11 12 13] [14 15 16 4] [17 18 19 5]};% [20 21 22 10] [2 23 24]};   % word PC 50 ignored bad trial
stimleg  = {["0", "25", "75"]  ["0", "25", "50", "75"] ["3" "5" "8" "100"] ["4" "6" "10" "100"]};
stimgrnames = {'Word Phase-coherence' 'Face Phase-coherence'    'Word Contrast'   'Face Contrast'};%   'Noise Con'   'Other'};

spect_base = nan(length(subjects),numchannels,200);

for zzz = 1:2
    for ccc =  elec_consider{zzz}
        S_temp = [];
        for ppp = 1:2
            for stimg = 1:length(stimgroups)
                for stimc = 1:length(stimgroups{stimg})
                    temp = mean(data(winlen_base, zzz, ccc, ppp, :, stimgroups{stimg}(stimc)), 5);
                    temp(isnan(temp)) = 0;
                    [Spectra,f_spectra] = getWaveletSpectrogram(temp, fsorig, [1, fupper], 'kjm');
                    S_temp = [S_temp; Spectra'];
                end
            end
        end
        spect_base (zzz,ccc,:) = mean(S_temp, 1);
    end
end

clear temp S_temp Spectra f_spectra

%% Plotting Raw Voltage + Spectra

stimgroups  = {5};    % only Face 100% Con/PC response
stimleg  = {"Face"};

% Window Length (0.5 s before and after stimulus) 
winlen_tc = find(epochtime >= -0.5 & epochtime <= 2.5);   % -0.5 s to 2.5 s


for zzz = 1:2
    for ccc = elec_consider{zzz}
        for ppp = 2     % only cat task
         
            sig_consider = squeeze(data(winlen_tc, zzz, ccc, ppp, :, stimgroups{1}));
            nTrials = size(sig_consider, 2);
            yMean = mean(sig_consider, 2);
            ySEM  = std(sig_consider, 0, 2) / sqrt(nTrials);
            ts = tinv([0.025, 0.975], nTrials - 1); % n-1 degrees of freedom -HH
            conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
%             if ppp == 1
%                 plot(smooth(yMean,40),'b-', 'LineWidth', 2);
%                 xline(find(epochtime == 0),'--k');
%                 fill([1:length(winlen_tc) length(winlen_tc):-1:1],conf, [0 0 .5],'EdgeColor',[0 0 .5],'FaceAlpha',.5,'HandleVisibility','off');
%             elseif ppp == 2
            subplot(4, 2, zzz)
            plot(smooth(yMean,40),'r-', 'LineWidth', 2);
            hold on;
            xline(find(epochtime == 0),'--k');
            fill([1:length(winlen_tc) length(winlen_tc):-1:1], smooth(conf,40), [.5 0 0],'EdgeColor',[.5 0 0],'FaceAlpha',.5,'HandleVisibility','off');
            xlim([winlen_tc(1) winlen_tc(end)]);
            ylabel('Normalized Amplitude');
            xlabel('Time relative to stimulus onset (s)');
            set(gca,'XTick',winlen_tc(1):1000:winlen_tc(end)); set(gca,'XTickLabel',epochtime(winlen_tc(1)):0.5:epochtime(winlen_tc(end)));
            %title(stimgrnames{stimg}+" "+stimleg{stimg}(stimc));     
            
            [Spectra,f_spectra] = getWaveletSpectrogram(yMean, fsorig, [1, fupper], 'kjm');
            spectra = log10(Spectra) - log10(squeeze(spect_base(zzz,ccc,:)));
            
            %plotting spectra
            figureprep([200 200 3800 2200]);
            subplot(4, 2, [2+zzz, 4+zzz, 6+zzz])
            load('dg_colormap.mat', 'cm') %KJM color scheme;
            uimagesc(winlen_tc,f_spectra,movmean(movmean(spectra, 1, 2), 1, 1),[-1.3 1.3]); axis xy; colormap(cm); colorbar
            xlim([winlen_tc(1) winlen_tc(end)]);
            ylabel('Frequency (Hz)');
            xlabel('Time (s)');
            set(gca,'XTick',winlen_tc(1):1000:winlen_tc(end)); set(gca,'XTickLabel',epochtime(winlen_tc(1)):0.5:epochtime(winlen_tc(end)));
            exportgraphics(gca, [sprintf( 'FvsC/eps/spectraSubj%d_ch%03d.eps', zzz, elec_consider{zzz} )]);
            %title(stimgrnames{stimg}+" "+stimleg{stimg}(stimc)); 
            
        end
    end
end

% Exporting only colorbar
% colormap(cm); colorbar; caxis([-1.3 1.3])
% exportgraphics(gca, [sprintf( 'FvsC/eps/spectraColorBar.eps' )]);


%% Broadband - Fixation v/s Categorization over con and phase coherence averaged over all visually responsive electrodes and face and word category (general overview)

gcc = {horzcat((9:11), 41, (44:46), (74:78)), horzcat((20:22), 28, 44, (86:87), 108, 110) };
pval = [];

stimgroups  = {[ 6      7      8      9       4  ] [10     11     12     13       5  ] [14    15     16      4   ] [17    18     19      5   ]};
stimleg     = {["0%", "25%", "50%", "75%", "100%"]                                     ["4%", "6%", "10%", "100%"]};
stimgrnames = {'Phase Coherence'                                                        'Contrast'};

winlen_tc = find(epochtime_bb >= -0.5 & epochtime_bb <= 2.5);   % -0.5 s to 2.5 s
winlen_sf = find(epochtime_bb >= 0 & epochtime_bb <= 1);        %    0 s to 1 s

fix_avgresp = nan(length(winlen_tc), length(gcc{1})+length(gcc{2}), 2*numreps, 2, 5); % fixation - time x number of visually responsive electrodes x 2*number of trials x # of categories X #of stim in each category
cat_avgresp = nan(length(winlen_tc), length(gcc{1})+length(gcc{2}), 2*numreps, 2, 5); % categorization

% Combining face and words categories over subjects and electrodes 
counter = 1;
for zzz = 1:2
    for ccc = gcc{zzz}
        for stimg = 1:2 % we need only 1st and 3rd group
            for stimc = 1 : length(stimgroups{2*stimg-1})
                fix_avgresp(:, counter,  1:6, stimg, stimc) = squeeze( bbdata_pc(winlen_tc, zzz, ccc, 1, :, stimgroups{2*stimg-1}(stimc))); % fixation 
                fix_avgresp(:, counter, 7:12, stimg, stimc) = squeeze( bbdata_pc(winlen_tc, zzz, ccc, 1, :, stimgroups{2*stimg}(stimc)));
                cat_avgresp(:, counter,  1:6, stimg, stimc) = squeeze( bbdata_pc(winlen_tc, zzz, ccc, 2, :, stimgroups{2*stimg-1}(stimc))); % categorization
                cat_avgresp(:, counter, 7:12, stimg, stimc) = squeeze( bbdata_pc(winlen_tc, zzz, ccc, 2, :, stimgroups{2*stimg}(stimc))); 
            end
        end
        counter = counter+1;
    end
end

assert(counter == length(gcc{1})+length(gcc{2})+1)

fix_avgresp = reshape(fix_avgresp, [length(winlen_tc), (length(gcc{1})+length(gcc{2}))*2*numreps, 2, 5]);
cat_avgresp = reshape(cat_avgresp, [length(winlen_tc), (length(gcc{1})+length(gcc{2}))*2*numreps, 2, 5]);

fix_avgresp = movmean(fix_avgresp, 40, 1);
cat_avgresp = movmean(cat_avgresp, 40, 1);

% Plot

for stimg = 1:2
    
    figureprep([200 200 3800 2200]);
    mx = -Inf; mn = Inf;
    
    for stimc = 1 : length(stimgroups{2*stimg-1})
        
        subplot(2,length(stimgroups{2*stimg-1}),stimc); hold on;
        
        nTrials = size(fix_avgresp(:,:,stimg,stimc),2);
        
        % Fixation
        yMean = mean(fix_avgresp(:,:,stimg,stimc),2);
        ySEM  = std(fix_avgresp(:,:,stimg,stimc), 0, 2)/(sqrt(nTrials));
        ts = tinv([0.025, 0.975], nTrials - 1);         % n-1 degrees of freedom -HH
        conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
        plot(yMean, 'b-', 'LineWidth', 2);
        xline(find(epochtime_bb == 0), '--k')
        fill([1:length(winlen_tc) length(winlen_tc):-1:1],conf, [0 0 .5],'EdgeColor','none','FaceAlpha',.3,'HandleVisibility','off');
        
        % Categorization
        yMean = mean(cat_avgresp(:,:,stimg,stimc),2);
        ySEM  = std(cat_avgresp(:,:,stimg,stimc), 0, 2)/(sqrt(nTrials));
        ts = tinv([0.025, 0.975], nTrials - 1);
        conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
        plot(yMean, 'r-', 'LineWidth', 2);
        fill([1:length(winlen_tc) length(winlen_tc):-1:1],conf, [.5 0 0],'EdgeColor','none','FaceAlpha',.3,'HandleVisibility','off');
        
        ax = axis;
        mx = max(mx,max(ax(3:4)));
        mn = min(mn,min(ax(3:4)));
        xlim([winlen_tc(1) winlen_tc(end)]);
        ylabel('Broadband Power (x-fold increase)');
        xlabel('Time relative to stimulus onset (s)');
        set(gca,'XTick',winlen_tc(1):100:winlen_tc(end)); set(gca,'XTickLabel',epochtime_bb(winlen_tc(1)):0.5:epochtime_bb(winlen_tc(end)));
        title(stimgrnames{stimg}+" "+stimleg{stimg}(stimc));
        legend('Fixation','Categorization');
        
    end
    
    for stimc = 1:length(stimgroups{2*stimg-1})
        
        subplot(2,length(stimgroups{2*stimg-1}),stimc);
        ax = axis;
        axis([ax(1:2) 0 5]);
        %exportgraphics(gca,[sprintf('FvsC/eps/STC%s%s.eps',stimgrnames{stimg},stimleg{stimg}(stimc))], 'BackgroundColor','none','ContentType','vector');
    
    end
        
    mx = -Inf; mn = Inf;
    
    for stimc = 1 : length(stimgroups{2*stimg-1})
        
        % Scaling Factor
        subplot(2,length(stimgroups{2*stimg-1}), length(stimgroups{2*stimg-1})+stimc); hold on;
        
        scaling_factor = zerodiv(squeeze(cat_avgresp(winlen_sf,:,stimg,stimc))+1 , squeeze(fix_avgresp(winlen_sf,:,stimg,stimc))+1, 1);     % zerodiv: knkutil % adding 1 so as now the broadband is x-fold and not increases
        
        nTrials = size(scaling_factor,2);
        yMean = mean(scaling_factor,2);
        ySEM  = std(scaling_factor,0,2)/(sqrt(nTrials));
        ts = tinv([0.025, 0.975], nTrials - 1);
        conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
        plot(yMean, 'k-', 'LineWidth', 2);
        yline(1, '--k')
        fill([1:length(winlen_sf) length(winlen_sf):-1:1],conf, [0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',.3,'HandleVisibility','off');
        
        for tt = 1: length(winlen_sf)    % CoCoSys Lab
            [~, pval(tt)] = ttest(squeeze(cat_avgresp(winlen_sf(tt),:,stimg,stimc)), squeeze(fix_avgresp(winlen_sf(tt),:,stimg,stimc)));
        end
        % convert to logical
        signific = nan(1, length(winlen_sf)); 
        signific(pval < 0.05) = 1;
        plot(1:length(winlen_sf), signific * 0.8, '.g', 'MarkerSize', 30);
        
        ax = axis;
        mx = max(mx,max(ax(3:4)));
        mn = min(mn,min(ax(3:4)));
        xlim([1 length(winlen_sf)]);
        ylabel('Scaling Factor');
        xlabel('Time relative to stimulus onset (s)');
        set(gca,'XTick',1:100:length(winlen_sf)); set(gca,'XTickLabel',epochtime_bb(winlen_sf(1)):0.5:epochtime_bb(winlen_sf(end)));
        
    end
    
    for stimc = 1:length(stimgroups{2*stimg-1})
        
        subplot(2,length(stimgroups{2*stimg-1}), length(stimgroups{2*stimg-1})+stimc);
        ax = axis;
        axis([ax(1:2) 0.6 1.8]);
        %exportgraphics(gca,[sprintf('FvsC/eps/SF%s%s.eps',stimgrnames{stimg},stimleg{stimg}(stimc))],'BackgroundColor','none','ContentType','vector');

    end
    figurewrite(sprintf('%s',stimgrnames{stimg}),[],[],'FvsC');
    %exportgraphics(gca,[sprintf('FvsC/eps/%s.eps',stimgrnames{stimg})]);
end


%% Alpha - Fixation v/s Categorization over con and phase coherence averaged over all visually responsive electrodes and face and word category (general overview)

gcc = {horzcat((9:11), 41, (44:46), (74:78)), horzcat((20:22), 28, 44, (86:87), 108, 110) };
pval = [];

stimgroups  = {[ 6      7      8      9       4  ] [10     11     12     13       5  ] [14    15     16      4   ] [17    18     19      5   ]};
stimleg     = {["0%", "25%", "50%", "75%", "100%"]                                     ["4%", "6%", "10%", "100%"]};
stimgrnames = {'Phase Coherence'                                                        'Contrast'};

winlen_tc = find(epochtime_bb >= -0.5 & epochtime_bb <= 2.5);   % -0.5 s to 2.5 s
winlen_sf = find(epochtime_bb >= 0 & epochtime_bb <= 2);        %    0 s to 1 s

fix_avgresp = zeros(length(winlen_tc), length(gcc{1})+length(gcc{2}), 2*numreps, 2, 5); % fixation - time x number of visually responsive electrodes x 2*number of trials x #of stim in each category
cat_avgresp = zeros(length(winlen_tc), length(gcc{1})+length(gcc{2}), 2*numreps, 2, 5); % categorization

% Combining face and words categories over subjects and electrodes 
counter = 1;
for zzz = 1:2
    for ccc = gcc{zzz}
        for stimg = 1:2 % we need only 1st and 3rd group
            for stimc = 1 : length(stimgroups{2*stimg-1})
                fix_avgresp(:, counter,  1:6, stimg, stimc) = squeeze( alphadata_pc(winlen_tc, zzz, ccc, 1, :, stimgroups{2*stimg-1}(stimc))); % fixation 
                fix_avgresp(:, counter, 7:12, stimg, stimc) = squeeze( alphadata_pc(winlen_tc, zzz, ccc, 1, :, stimgroups{2*stimg}(stimc)));
                cat_avgresp(:, counter,  1:6, stimg, stimc) = squeeze( alphadata_pc(winlen_tc, zzz, ccc, 2, :, stimgroups{2*stimg-1}(stimc))); % categorization
                cat_avgresp(:, counter, 7:12, stimg, stimc) = squeeze( alphadata_pc(winlen_tc, zzz, ccc, 2, :, stimgroups{2*stimg}(stimc))); 
            end
        end
        counter = counter+1;
    end
end

assert(counter == length(gcc{1})+length(gcc{2})+1)

fix_avgresp = reshape(fix_avgresp, [length(winlen_tc), (length(gcc{1})+length(gcc{2}))*2*numreps, 2, 5]);
cat_avgresp = reshape(cat_avgresp, [length(winlen_tc), (length(gcc{1})+length(gcc{2}))*2*numreps, 2, 5]);

fix_avgresp = movmean(fix_avgresp, 40, 1);
cat_avgresp = movmean(cat_avgresp, 40, 1);

% Plot

for stimg = 1:2
    
    figureprep([200 200 3800 2200]);
    mx = -Inf; mn = Inf;
    
    for stimc = 1 : length(stimgroups{2*stimg-1})
        
        subplot(2,length(stimgroups{2*stimg-1}),stimc); hold on;
        
        nTrials = size(fix_avgresp(:,:,stimg,stimc),2);
        
        % Fixation
        yMean = mean(fix_avgresp(:,:,stimg,stimc),2);
        ySEM  = std(fix_avgresp(:,:,stimg,stimc), 0, 2)/(sqrt(nTrials));
        ts = tinv([0.025, 0.975], nTrials - 1);         % n-1 degrees of freedom -HH
        conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
        plot(yMean, 'b-', 'LineWidth', 2);
        xline(find(epochtime_bb == 0), '--k')
        yline(0, '--k');
        fill([1:length(winlen_tc) length(winlen_tc):-1:1],conf, [0 0 .5],'EdgeColor','none','FaceAlpha',.3,'HandleVisibility','off');
        
        % Categorization
        yMean = mean(cat_avgresp(:,:,stimg,stimc),2);
        ySEM  = std(cat_avgresp(:,:,stimg,stimc), 0, 2)/(sqrt(nTrials));
        ts = tinv([0.025, 0.975], nTrials - 1);
        conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
        plot(yMean, 'r-', 'LineWidth', 2);
        fill([1:length(winlen_tc) length(winlen_tc):-1:1],conf, [.5 0 0],'EdgeColor','none','FaceAlpha',.3,'HandleVisibility','off');
        
        ax = axis;
        mx = max(mx,max(ax(3:4)));
        mn = min(mn,min(ax(3:4)));
        xlim([winlen_tc(1) winlen_tc(end)]);
        ylabel('Broadband Power (x-fold increase)');
        xlabel('Time relative to stimulus onset (s)');
        set(gca,'XTick',winlen_tc(1):100:winlen_tc(end)); set(gca,'XTickLabel',epochtime_bb(winlen_tc(1)):0.5:epochtime_bb(winlen_tc(end)));
        title(stimgrnames{stimg}+" "+stimleg{stimg}(stimc));
        legend('Fixation','Categorization');
        
    end
    
    for stimc = 1:length(stimgroups{2*stimg-1})
        
        subplot(2,length(stimgroups{2*stimg-1}),stimc);
        ax = axis;
        axis([ax(1:2) -1 1]);
        %exportgraphics(gca, [sprintf('alphaFvsC/eps/STC%s%s.eps', stimgrnames{stimg}, stimleg{stimg}(stimc))], 'BackgroundColor', 'none', 'ContentType', 'vector');
    end
        
    mx = -Inf; mn = Inf;
    
    for stimc = 1 : length(stimgroups{2*stimg-1})
        
        % Scaling Factor
        subplot(2,length(stimgroups{2*stimg-1}), length(stimgroups{2*stimg-1})+stimc); hold on;
        
        scaling_factor = zerodiv(squeeze(cat_avgresp(winlen_sf,:,stimg,stimc))+1 , squeeze(fix_avgresp(winlen_sf,:,stimg,stimc))+1, 1);     % zerodiv: knkutil % adding 1 so as now the broadband is x-fold and not increases
        
        nTrials = size(scaling_factor,2);
        yMean = median(scaling_factor,2);
        ySEM  = std(scaling_factor,0,2)/(sqrt(nTrials));
        ts = tinv([0.025, 0.975], nTrials - 1);
        conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
        plot(yMean, 'k-', 'LineWidth', 2);
        yline(1, '--k')
        fill([1:length(winlen_sf) length(winlen_sf):-1:1],conf, [0.5 0.5 0.5],'EdgeColor','none','FaceAlpha',.3,'HandleVisibility','off');
        
        for tt = 1: length(winlen_sf)    % CoCoSys Lab
            [~, pval(tt)] = ttest(squeeze(cat_avgresp(winlen_sf(tt),:,stimg,stimc)), squeeze(fix_avgresp(winlen_sf(tt),:,stimg,stimc)));
        end
        % convert to logical
        signific = nan(1, length(winlen_sf)); 
        signific(pval < 0.05) = 1;
        plot(1:length(winlen_sf), signific * -0.2, '.g', 'MarkerSize', 30);
        
        ax = axis;
        mx = max(mx,max(ax(3:4)));
        mn = min(mn,min(ax(3:4)));
        xlim([1 length(winlen_sf)]);
        ylabel('Scaling Factor');
        xlabel('Time relative to stimulus onset (s)'); 
        set(gca,'XTick',1:100:length(winlen_sf)); set(gca,'XTickLabel',epochtime_bb(winlen_sf(1)):0.5:epochtime_bb(winlen_sf(end)));
        
    end
    
    for stimc = 1:length(stimgroups{2*stimg-1})
        
        subplot(2,length(stimgroups{2*stimg-1}), length(stimgroups{2*stimg-1})+stimc);
        ax = axis;
        axis([ax(1:2) -0.3 1.2]);
        %exportgraphics(gca, [sprintf('alphaFvsC/eps/SF%s%s.eps', stimgrnames{stimg}, stimleg{stimg}(stimc))], 'BackgroundColor', 'none', 'ContentType', 'vector');
    end
    figurewrite(sprintf('%s',stimgrnames{stimg}),[],[],'alphaFvsC');
    %exportgraphics(gca,[sprintf('FvsC/eps/%s.eps',stimgrnames{stimg})]);
end
%% Fixation v/s Categorization over con and phase coherence for electrode categires as well as face and word category (supp to general overview)

%  gcc = {horzcat((9:11), 15, 41, (44:46), (74:78)), horzcat((20:22), 28, 44, (86:87), 108, 110) };
%  plotname = 'FvsC_allVisGood';
%  gcc = {horzcat((9:11), 15, 41, (44:46)),          horzcat((20:22), 28, (86:87))};   % All-VTC
%  plotname = 'FvsC_allVTC';
%  gcc = {[74, 75], []};   %EVC
%  plotname = 'FvsC_EVC';
%  gcc = {[76, 77, 78], [44, 108, 110]};   % LVC
%  plotname = 'FvsC_LVC';
%  gcc =  {[9 15 44],[21 87]};      % VTC
%  plotname = 'FvsC_VTC';
%  gcc =  {[10 41 45 46],[28]};     % FFA
%  plotname = 'FvsC_FFA';
 gcc =  {[11], [20 22 86]};        % VWFA
 plotname = 'FvsC_VWFA';

stimgroups  = {[ 6      7      8      9       4  ] [10     11     12     13       5  ] [14    15     16      4   ] [17    18     19      5   ]};
stimleg     = {["0%", "25%", "50%", "75%", "100%"] ["0%", "25%", "50%", "75%", "100%"] ["3%", "5%",  "8%", "100%"] ["4%", "6%", "10%", "100%"]};
stimgrnames = {'Word Phase Coherence'               'Face Phase Coherence'               'Word Contrast'              'Face Contrast'};

winlen_tc = find(epochtime_bb >= -0.5 & epochtime_bb <= 2.5);   % -0.5 s to 2.5 s
winlen_sf = find(epochtime_bb >= 0 & epochtime_bb <= 1);        %    0 s to 1 s
winlen_pt = winlen_sf;      % win len for peak time

avgresp = nan(length(winlen_tc), length(gcc{1})+length(gcc{2}), numtasks, numreps, length(stimgroups), 5); %  time x number of visually responsive electrodes x tasks x trials x categories x stim
peaktime = nan(numtasks+2, length(stimgroups), 5);

% Extarcting only electrodes under consideration over both subjects
counter = 1;
for zzz = 1 : 2
    for ccc = gcc{zzz}
        for stimg = 1 : length(stimgroups)
            for stimc = 1 : length(stimgroups{stimg})
                avgresp(:, counter,  :, :, stimg, stimc) = squeeze(bbdata_pc(winlen_tc, zzz, ccc, :, :, stimgroups{stimg}(stimc)));
            end
        end
        counter = counter+1;
    end
end
 
assert(counter-1 == length(gcc{1})+length(gcc{2}));

avgresp = permute(avgresp, [1 2 4 3 5 6]);
avgresp = reshape(avgresp, [length(winlen_tc), (length(gcc{1})+length(gcc{2})) * numreps, numtasks, length(stimgroups), 5]);
avgresp = movmean(avgresp, 40, 1);

% Plot

for stimg = 1 : length(stimgroups)
    
    figureprep([200 200 3800 2200]);
    mx = -Inf; mn = Inf;
    
    for stimc = 1 : length(stimgroups{stimg})
        
        subplot(2,length(stimgroups{stimg}),stimc); hold on;
        
        nTrials = size(avgresp, 2);
        
        % Fixation
        yMean = mean(avgresp(:, :, 1, stimg, stimc), 2);
        [~, peaktime(1, stimg, stimc)] = max(yMean(winlen_pt));
        peaktime(1, stimg, stimc) = peaktime(1, stimg, stimc) + winlen_pt(1);
        ySEM  = std(avgresp(:, :, 1, stimg, stimc), 0, 2) / sqrt(nTrials);
        ts = tinv([0.025, 0.975], nTrials - 1); % n-1 degrees of freedom -HH
        conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
        plot(yMean, 'b-', 'LineWidth', 2);
        fill([1:length(winlen_tc) length(winlen_tc):-1:1],conf, [0 0 .5],'EdgeColor',[0 0 .5],'FaceAlpha',.5,'HandleVisibility','off');
        scatter(peaktime(1, stimg, stimc), yMean(peaktime(1, stimg, stimc)), 'g*');
        %errorbar(1:length(epochtime_bb),yMean,2*ySEM','k.');
        
        % Categorization
        yMean = mean(avgresp(:, :, 2, stimg, stimc), 2);
        [~, peaktime(2, stimg, stimc)] = max(yMean(winlen_pt));
        peaktime(2, stimg, stimc) = peaktime(2, stimg, stimc) + winlen_pt(1);
        ySEM  = std(avgresp(:, :, 2, stimg, stimc), 0, 2) / sqrt(nTrials);
        ts = tinv([0.025, 0.975], nTrials - 1); % n-1 degrees of freedom -HH
        conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
        scatter(peaktime(2, stimg, stimc), yMean(peaktime(2, stimg, stimc)), 'g*');
        plot(yMean, 'r-', 'LineWidth', 2);
        fill([1:length(winlen_tc) length(winlen_tc):-1:1],conf, [.5 0 0],'EdgeColor',[.5 0 0],'FaceAlpha',.5,'HandleVisibility','off');
        scatter(peaktime(2, stimg, stimc), yMean(peaktime(2, stimg, stimc)), 'g*');
        %errorbar(1:length(epochtime_bb),yMean,2*ySEM','k.');
       
        ax = axis;
        mx = max(mx,max(ax(3:4)));
        mn = min(mn,min(ax(3:4)));
        xlim([1 length(winlen_tc)]);
        ylabel('BB response (x-fold)');
        xlabel('t (s)');
        set(gca,'XTick',1:100:length(winlen_tc)); set(gca,'XTickLabel',epochtime_bb(winlen_tc(1)):0.5:epochtime_bb(winlen_tc(end)));
        title(stimgrnames{stimg}+" "+stimleg{stimg}(stimc));
        legend('F','C');
        
    end
    
    for stimc = 1:length(stimgroups{stimg})
        
        subplot(2,length(stimgroups{stimg}),stimc); hold on;
        ax = axis;
        axis([ax(1:2) mn mx]);
        %exportgraphics(gca,[sprintf('FvsC/%s%s/',stimgrnames{stimg},stimleg{stimg}(stimc)) 'STC.eps']);
    
    end
        
    mx = -Inf; mn = Inf;
    
    for stimc = 1 : length(stimgroups{stimg})
        
        subplot(2,length(stimgroups{stimg}), length(stimgroups{stimg})+stimc); hold on;
        
        scaling_factor = squeeze(avgresp(winlen_sf, :, 2, stimg, stimc)+1)./ squeeze(avgresp(winlen_sf, :, 1, stimg, stimc)+1);
        % scaling_factor = scaling_factor./(squeeze(fix_avgresp(:,:,stimg,stimc)) - 1);
        
        
        nTrials = size(scaling_factor, 2);
        yMean = mean(scaling_factor, 2);
        abvThr = ge(yMean, 1.2);
        i_abvThr = find(diff([0 abvThr' 0]) == 1);
        temp = i_abvThr(find(diff([0 abvThr' 0]) == -1) - i_abvThr >= 20);
        %peaktime(4, stimg, stimc) = temp(1);
        [~, peaktime(3, stimg, stimc)] = max(yMean);
        ySEM  = std(scaling_factor, 0, 2) / sqrt(nTrials);
        ts = tinv([0.025, 0.975], nTrials - 1); % n-1 degrees of freedom -HH
        conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
        plot(yMean,'k-', 'LineWidth', 2);
        fill([1:length(winlen_sf) length(winlen_sf):-1:1],conf, [.5 0.5 0.5],'EdgeColor',[.5 0.5 0.5],'FaceAlpha',.5,'HandleVisibility','off');
        scatter(peaktime(3, stimg, stimc), yMean(peaktime(3, stimg, stimc)), 'g*');
        %scatter(peaktime(4, stimg, stimc), yMean(peaktime(4, stimg, stimc)), 'c*');
        %[h,p,ci,stats] = ttest2(mean(stimresp_c(:,:,2,stimg,stimc),2),mean(avgresp(:,:,2,stimg,stimc),2));
%         for t = 1:length(yMean)     % CoCoSys Lab
%             [~, pval(t)] = ttest(squeeze(cat_avgresp(t,:,stimg,stimc)), squeeze(fix_avgresp(t,:,stimg,stimc)));
%         end
%         % convert to logical
%         signific = nan(1, length(yMean)); 
%         signific(pval < 0.05) = 1;
%         plot(signific * -0.03, '.k');
        % indicate what we're showing
        %text(10.2, -0.03, 'p < 0.001');
        ax = axis;
        axis([100 300 ax(3:4)]);
        mx = max(mx,max(ax(3:4)));
        mn = min(mn,min(ax(3:4)));
        xlim([1 length(winlen_sf)]);
        ylabel('BB response (x-fold)');
        xlabel('t (s)');
        set(gca,'XTick',1:100:length(winlen_sf)); set(gca,'XTickLabel',epochtime_bb(winlen_sf(1)):0.5:epochtime_bb(winlen_sf(end)));
        %title(stimgrnames{stimg}+" "+stimleg{stimg}(stimc));
        title('Scaling Factor');
        
    end
    
    for stimc = 1:length(stimgroups{stimg})
        
        subplot(2,length(stimgroups{stimg}), length(stimgroups{stimg})+stimc); hold on;
        ax = axis;
        axis([ax(1:2) mn mx]);
        %exportgraphics(gca,[sprintf('FvsC/%s%s/',stimgrnames{stimg},stimleg{stimg}(stimc)) 'SF.eps']);

    end
    figurewrite(sprintf('%s',stimgrnames{stimg}),[],[],plotname);
end

peaktime_VWFA = peaktime;

%% Peak Times

% peaktime_all = peaktime;
% peaktime_EVC = peaktime;
% peaktime_LVC = peaktime;
% peaktime_VTC = peaktime;
% peaktime_FFA = peaktime;
% peaktime_VWFA = peaktime;
% 
% for stimg = 1 : length(stimgroups)
% 
%         subplot(2,length(stimgroups), stimg); hold on;
%         plot(squeeze(peaktime_all(2, stimg, :) - peaktime(1, stimg, :)));
%         plot(squeeze(peaktime_EVC(2, stimg, :) - peaktime(1, stimg, :)));
%         plot(squeeze(peaktime_LVC(2, stimg, :) - peaktime(1, stimg, :)));
%         plot(squeeze(peaktime_VTC(2, stimg, :) - peaktime(1, stimg, :)));
%         plot(squeeze(peaktime_FFA(2, stimg, :) - peaktime(1, stimg, :)));
%         plot(squeeze(peaktime_VWFA(2, stimg, :) - peaktime(1, stimg, :)));
% 
%         legend('AllVisGood', 'EVC', 'LVC', 'VTC', 'FFA', 'VWFA')
% 
% end

figure
% for ttt = 1:numtasks
%     for stimg = 1 : length(stimgroups)
% 
%             subplot(2,length(stimgroups), (ttt-1) * length(stimgroups) + stimg); hold on;
%             bar(squeeze(peaktime_all(ttt, stimg, :) ));
%             bar(squeeze(peaktime_EVC(ttt, stimg, :) ));
%             bar(squeeze(peaktime_LVC(ttt, stimg, :) ));
%             bar(squeeze(peaktime_VTC(ttt, stimg, :) ));
%             bar(squeeze(peaktime_FFA(ttt, stimg, :) ));
%             bar(squeeze(peaktime_VWFA(ttt, stimg, :)));
% 
%             legend('AllVisGood', 'EVC', 'LVC', 'VTC', 'FFA', 'VWFA')
% 
%     end
% end
% 
% bar([squeeze(peaktime_all(ttt, stimg, :) ); squeeze(peaktime_EVC(ttt, stimg, :) ); squeeze(peaktime_LVC(ttt, stimg, :) ); squeeze(peaktime_VTC(ttt, stimg, :) ); squeeze(peaktime_FFA(ttt, stimg, :) ); squeeze(peaktime_VWFA(ttt, stimg, :))]);

for ttt = 1:numtasks
    for stimg = 1 : length(stimgroups)

            subplot(2,length(stimgroups), (ttt-1) * length(stimgroups) + stimg); hold on;
            plot(squeeze(peaktime_all(ttt, stimg, :) ));
            plot(squeeze(peaktime_EVC(ttt, stimg, :) ));
            plot(squeeze(peaktime_LVC(ttt, stimg, :) ));
            plot(squeeze(peaktime_VTC(ttt, stimg, :) ));
            plot(squeeze(peaktime_FFA(ttt, stimg, :) ));
            plot(squeeze(peaktime_VWFA(ttt, stimg, :)));

            legend('AllVisGood', 'EVC', 'LVC', 'VTC', 'FFA', 'VWFA')

    end
end

        
%% BB_PC averaged over time window - Different stim, categoreies, and electrode locations

 gcc = {horzcat((9:11), 15, 41, (44:46), (74:78)), horzcat((20:22), 28, 44, (86:87), 108, 110) }; % All Visually Good
 plotname = 'allVisGood';
%  gcc = {horzcat((9:11), 15, 41, (44:46)),          horzcat((20:22), 28, (86:87))};   % All-VTC
%  plotname = 'allVTC';
%  gcc = {[74, 75], []};                   % EVC
%  plotname = 'EVC';
%  gcc = {[76, 77, 78], [44, 108, 110]};   % LVC
%  plotname = 'LVC';
%  gcc =  {[9 15 44],[21 87]};             % VTC
%  plotname = 'VTC';
%  gcc =  {[10 41 45 46],[28]};            % FFA
%  plotname = 'FFA';
%  gcc =  {[11], [20 22 86]};              % VWFA
%  plotname = 'VWFA';

stimgroups  = {[ 6      7      8      9       4  ] [10     11     12     13       5  ] [14    15     16      4   ] [17    18     19      5   ]};
stimleg     = {["0%", "25%", "50%", "75%", "100%"] ["0%", "25%", "50%", "75%", "100%"] ["3%", "5%",  "8%", "100%"] ["4%", "6%", "10%", "100%"]};
stimgrnames = {'Word Phase Coherence'               'Face Phase Coherence'               'Word Contrast'              'Face Contrast'};

% Time window for estimating responsiveness 100 - 500 ms
winlen_resp = find(epochtime_bb > 0.1 & epochtime_bb < 0.5);

avgresp = zeros(length(gcc{1})+length(gcc{2}), numreps, numtasks, length(stimgroups), 5); % number of visually responsive electrodes per region x number of trials x # of tasks (fix/cat) x # of categories x stim in each category

% Plot

xxx = [];
yyy = [];

counter = 1;
for zzz = 1:2
    for ccc = gcc{zzz}
        for stimg = 1:length(stimgroups)
            for stimc = 1 : length(stimgroups{stimg})
                avgresp(counter, :, :, stimg, stimc) = squeeze(mean(bbdata_pc(winlen_resp, zzz, ccc, :, :, stimgroups{stimg}(stimc)), 1))'; 
            end
        end
        counter = counter+1;
    end
end

avgresp = reshape(avgresp, [(length(gcc{1})+length(gcc{2})) * numreps, numtasks, length(stimgroups), 5]);

posx = 1;
figure
mx = -Inf; mn = Inf;

    
for stimg = 1:length(stimgroups)

    temp_yMean =[];
    temp_posx  =[];
    
    for stimc = 1 : length(stimgroups{stimg})
        
        subplot(3, 2, [3, 5]); hold on;
        
        nTrials = size(avgresp(:, :, stimg, stimc), 1);
        
        % Fixation
        yMean = mean(avgresp(:, 1, stimg, stimc), 1);
        ySEM  =  std(avgresp(:, 1, stimg, stimc), 0, 1) / sqrt(nTrials);
        ts = tinv([0.025, 0.975], nTrials - 1); % n-1 degrees of freedom -HH
        conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
        bar(posx, yMean, 'FaceColor', [0 0 0.5], 'FaceAlpha', .5, EdgeColor='none');
        errorbar(posx, yMean, 2*ySEM','k.');
        
        ax = axis;
        mx = max(mx,max(ax(3:4)));
        mn = min(mn,min(ax(3:4)));
        title('Fixation');
        xticks([])
        exportgraphics(gca,[sprintf('FvsC/Avg/%s/Fix.eps', plotname)], 'BackgroundColor', 'none', 'ContentType', 'vector');
 
        
        subplot(3, 2, [4, 6]); hold on;
        
        % Categorization
        yMean = mean(avgresp(:, 2, stimg, stimc), 1);
        ySEM  =  std(avgresp(:, 2, stimg, stimc), 0, 1) / sqrt(nTrials);
        ts = tinv([0.025, 0.975], nTrials - 1); % n-1 degrees of freedom -HH
        conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
        bar(posx, yMean, 'FaceColor', [0.5 0 0], 'FaceAlpha', .5, EdgeColor='none');
        errorbar(posx, yMean, 2*ySEM','k.');
        
        ax = axis;
        mx = max(mx,max(ax(3:4)));
        mn = min(mn,min(ax(3:4)));
        title('Categorization');
        xticks([])
        exportgraphics(gca,[sprintf('FvsC/Avg/%s/Cat.eps', plotname)], 'BackgroundColor', 'none', 'ContentType', 'vector');
        
        
        subplot(3, 2, 2); hold on;
        
        scaling_factor = (squeeze(avgresp(:, 2, stimg, stimc)) + 1) ./ (squeeze(avgresp(:, 1, stimg, stimc)) + 1);
        %scaling_factor = scaling_factor./(squeeze(fix_avgresp(:,:,stimg,stimc)) - 1);
        
        nTrials = size(scaling_factor, 1);
        yMean = median(scaling_factor);
        temp_yMean = [temp_yMean, yMean];
        temp_posx = [temp_posx posx];
        ySEM  = std(scaling_factor) / (sqrt(nTrials));
        errorbar(posx, yMean, 2*ySEM', 'k.');
        

       
%         yMean1 = 0.5 * (nanmean(reactiontime(2,:,stimgroups{stimg}(stimc)),2) + nanmean(reactiontime(1,:,stimgroups{stimg}(stimc)),2));
%         xxx = [xxx yMean];
%         yyy = [yyy yMean1];
       
        % errorbar(posx, yMean, 2*ySEM','k.');
        [~, pval] = ttest(squeeze(avgresp(:, 2, stimg, stimc)), squeeze(avgresp(:, 1, stimg, stimc)));
        % convert to logical
        signific = nan; 
        signific(pval < 0.05) = 1;
        plot(posx, signific * 2, '.g');

        
        posx = posx + 1;
        
    end
    
    plot(temp_posx, temp_yMean, 'c', LineWidth = 2);
    yline(1, '--k');
    xticks([])
    exportgraphics(gca,[sprintf('FvsC/Avg/%s/SF.eps', plotname)], 'BackgroundColor', 'none', 'ContentType', 'vector');

    %plot reaction times
%     subplot(3, 2, [5 , 6]);
%     hold on;



%     for stimc = 1:length(stimgroups{stimg})
%         yMean = 0.5 * (nanmean(reactiontime(2,:,stimgroups{stimg}(stimc)),2) + nanmean(reactiontime(1,:,stimgroups{stimg}(stimc)),2));
%         y = [y, yMean];
%         bar(posx,yMean);
%         posx = posx + 1;
%     end
    
    posx = posx + 2;

    
end

clear temp_posx temp_yMean


% 
% figure, scatter(xxx, yyy)
% corr(xxx', yyy')

%%
x1=[];
x2=[];
y1=[];
y2=[];
x= [];
y =[];
posx = 1;
figure
mx = -Inf; mn = Inf;
    
for stimg = 1:length(stimgroups)
    
    for stimc = 1 : length(stimgroups{stimg})
        
        subplot(3, 2, 1); hold on;
        
        nTrials = size(avgresp(:, :, stimg, stimc), 1);
        
        % Fixation
        yMean = mean(avgresp_vtc(:, 1, stimg, stimc), 1) - mean(avgresp_lvc(:, 1, stimg, stimc), 1) - mean(avgresp_evc(:, 1, stimg, stimc), 1);
        x1 = [x1, mean(avgresp_lvc(:, 1, stimg, stimc), 1)];
        y1 = [y1, mean(avgresp_evc(:, 1, stimg, stimc), 1)];
        ySEM  =  std(avgresp_vtc(:, 1, stimg, stimc), 0, 1) / sqrt(nTrials) - std(avgresp_lvc(:, 1, stimg, stimc), 0, 1) / sqrt(nTrials) - std(avgresp_evc(:, 1, stimg, stimc), 0, 1) / sqrt(nTrials);
        ts = tinv([0.025, 0.975], nTrials - 1); % n-1 degrees of freedom -HH
        conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
        bar(posx, yMean,'w');
        errorbar(posx, yMean, 2*ySEM','k.');
        %fill([1:length(epochtime_bb) length(epochtime_bb):-1:1],conf, [0 0 .5],'EdgeColor',[0 0 .5],'FaceAlpha',.5,'HandleVisibility','off');
        ax = axis;
        mx = max(mx,max(ax(3:4)));
        mn = min(mn,min(ax(3:4)));
        title('Fixation');
        
        
        subplot(3, 2, 2); hold on;
        
        % Categorization
        yMean = mean(avgresp_vtc(:, 2, stimg, stimc), 1) - mean(avgresp_lvc(:, 2, stimg, stimc), 1) - mean(avgresp_evc(:, 2, stimg, stimc), 1);
        x2 = [x2, mean(avgresp_lvc(:, 2, stimg, stimc), 1)];
        y2 = [y2, mean(avgresp_evc(:, 2, stimg, stimc), 1)];
        ySEM  =  std(avgresp_vtc(:, 2, stimg, stimc), 0, 1) / sqrt(nTrials) - std(avgresp_lvc(:, 2, stimg, stimc), 0, 1) / sqrt(nTrials) - std(avgresp_evc(:, 2, stimg, stimc), 0, 1) / sqrt(nTrials);
        ts = tinv([0.025, 0.975], nTrials - 1); % n-1 degrees of freedom -HH
        conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
        bar(posx, yMean,'w');
        errorbar(posx, yMean, 2*ySEM','k.');
        %fill([1:length(epochtime_bb) length(epochtime_bb):-1:1],conf, [0 0 .5],'EdgeColor',[0 0 .5],'FaceAlpha',.5,'HandleVisibility','off');
        ax = axis;
        mx = max(mx,max(ax(3:4)));
        mn = min(mn,min(ax(3:4)));
        title('Categorization');
        
        subplot(3, 2, [3 , 4]); hold on;
        
        scaling_factor = squeeze(mean(avgresp_vtc(:, 2, stimg, stimc), 1) - mean(avgresp_lvc(:, 2, stimg, stimc), 1) - mean(avgresp_evc(:, 2, stimg, stimc), 1)) ./ squeeze(mean(avgresp_vtc(:, 1, stimg, stimc), 1) - mean(avgresp_lvc(:, 1, stimg, stimc), 1) - mean(avgresp_evc(:, 1, stimg, stimc), 1));
        %scaling_factor = scaling_factor./(squeeze(fix_avgresp(:,:,stimg,stimc)) - 1);
        
        nTrials = size(scaling_factor, 1);
        yMean = mean(scaling_factor);
        x = [x, yMean];
        ySEM  = std(scaling_factor) / (sqrt(nTrials));
        bar(posx, yMean,'w');
        plot(posx, yMean, '*g');
        %errorbar(posx, yMean, 2*ySEM','k.');
        [~, pval] = ttest(squeeze(avgresp(:, 2, stimg, stimc)), squeeze(avgresp(:, 1, stimg, stimc)));
        % convert to logical
        signific = nan; 
        signific(pval < 0.05) = 1;
        plot(posx, signific * -0.03, '.k');

        
        posx = posx + 1;
        
    end
    

    %plot reaction times
    subplot(3, 2, [5 , 6]);
    hold on;

    for stimc = 1:length(stimgroups{stimg})
        yMean = 0.5 * (nanmean(reactiontime(2,:,stimgroups{stimg}(stimc)),2) + nanmean(reactiontime(1,:,stimgroups{stimg}(stimc)),2));
        y = [y, yMean];
        bar(posx,yMean);
        posx = posx + 1;
    end
    posx = posx + 1;

    
end

%%  evc vtc lvc sf 



% gcc = {horzcat((9:11), 15, 41, (44:46), (74:78)), horzcat((20:22), 28, 44, (86:87), 108, 110) };
%  gcc = {[74, 75], []};   %EVC
% gcc = {[76, 77, 78], [44, 108, 110]};   % LVC
% gcc =  {[9 15 44],[21 87]};   % VTCns
% gcc =  {[10 41 45 46],[28]};    % FFA
% gcc =  {[11], [20 22 86]} % VWFA

% gcc = {horzcat((9:11), 15, 41, (44:46)), horzcat((20:22), 28,(86:87))}; %VTC

stimgroups  = {[ 6      7      8      9       4  ] [10     11     12     13       5  ] [14    15     16      4   ] [17    18     19      5   ]};
stimleg     = {["0%", "25%", "50%", "75%", "100%"] ["0%", "25%", "50%", "75%", "100%"] ["4%", "6%", "10%", "100%"] ["4%", "6%", "10%", "100%"]};
stimgrnames = {'Word Phase Coherence'               'Face Phase Coherence'               'Word Contrast'              'Face Contrast'};

winlen = 1:601; % -0.5 to 2.5 s

avgresp = zeros(length(epochtime_bb), length(gcc{1})+length(gcc{2}), numtasks, numreps, length(stimgroups), 5); % fixation - time x number of visually responsive electrodes x tasks x trials x categories x stim
peaktime = zeros(numtasks, length(stimgroups), 5);

% Extarcting only good electrodes over both subjects
counter = 1;
for zzz = 1 : 2
    for ccc = gcc{zzz}
        for stimg = 1 : length(stimgroups)
            for stimc = 1 : length(stimgroups{stimg})
                avgresp(:, counter,  :, :, stimg, stimc) = squeeze(bb_data(:, zzz, ccc, :, :, stimgroups{stimg}(stimc)));
            end
        end
        counter = counter+1;
    end
end
 
assert(counter-1 == length(gcc{1})+length(gcc{2}));

avgresp = permute(avgresp, [1 2 4 3 5 6]);
avgresp = reshape(avgresp, [length(epochtime_bb), (length(gcc{1})+length(gcc{2})) * numreps, numtasks, length(stimgroups), 5]);

% Plot
%%
for stimg = 1 : length(stimgroups)
    
    figure
    mx = -Inf; mn = Inf;
    
    for stimc = 1 : length(stimgroups{stimg})
        
        subplot(2,length(stimgroups{stimg}),stimc); hold on;
        
        nTrials = size(avgresp, 2);
        
        % Fixation
        yMean = smooth(mean(avgresp_lvc(:, :, 2, stimg, stimc), 2), 40);
        [~, peaktime(1, stimg, stimc)] = max(yMean(find(epochtime_bb > 0 & epochtime_bb < 0.5)));
        ySEM  = smooth(std(avgresp(:, :, 1, stimg, stimc), 0, 2) / sqrt(36*nTrials), 40);
        ts = tinv([0.025, 0.975], nTrials - 1); % n-1 degrees of freedom -HH
        conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
        plot(yMean, 'b-', 'LineWidth', 2);
        %fill([1:length(epochtime_bb) length(epochtime_bb):-1:1],conf, [0 0 .5],'EdgeColor',[0 0 .5],'FaceAlpha',.5,'HandleVisibility','off');
        %errorbar(1:length(epochtime_bb),yMean,2*ySEM','k.');
        
        % Categorization
        yMean = smooth(mean(avgresp_vtc(:, :, 2, stimg, stimc), 2), 40);
        [~, peaktime(2, stimg, stimc)] = max(yMean(find(epochtime_bb > 0 & epochtime_bb < 0.5)));
        ySEM  = smooth(std(avgresp(:, :, 2, stimg, stimc), 0, 2) / sqrt(36*nTrials), 40);
        ts = tinv([0.025, 0.975], nTrials - 1); % n-1 degrees of freedom -HH
        conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
        plot(yMean, 'r-', 'LineWidth', 2);
        %fill([1:length(epochtime_bb) length(epochtime_bb):-1:1],conf, [.5 0 0],'EdgeColor',[.5 0 0],'FaceAlpha',.5,'HandleVisibility','off');
        %errorbar(1:length(epochtime_bb),yMean,2*ySEM','k.');
       
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
    
    for stimc = 1:length(stimgroups{stimg})
        
        subplot(2,length(stimgroups{stimg}),stimc); hold on;
        ax = axis;
        axis([ax(1:2) mn mx]);
        %exportgraphics(gca,[sprintf('FvsC/%s%s/',stimgrnames{stimg},stimleg{stimg}(stimc)) 'STC.eps']);
    
    end
        
    mx = -Inf; mn = Inf;
    
    for stimc = 1 : length(stimgroups{stimg})
        
        subplot(2,length(stimgroups{stimg}), length(stimgroups{stimg})+stimc); hold on;
        
        scaling_factor = zerodiv(squeeze(mean(avgresp_lvc(:,:,2,stimg,stimc),2)) , squeeze(mean(avgresp_vtc(:,:,2,stimg,stimc),2)), 1);%...
%             - zerodiv(squeeze(mean(avgresp_lvc(:,:,1,stimg,stimc),2)) , squeeze(mean(avgresp_vtc(:,:,1,stimg,stimc),2)), 1)
        %scaling_factor = scaling_factor./(squeeze(fix_avgresp(:,:,stimg,stimc)) - 1);
        
        nTrials = size(scaling_factor, 2);
        yMean = smooth(scaling_factor,40);
        ySEM  = smooth(std(scaling_factor, 0, 2) / sqrt(36*nTrials), 40);
        ts = tinv([0.025, 0.975], nTrials - 1); % n-1 degrees of freedom -HH
        conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
        plot(yMean,'k-', 'LineWidth', 2);
        %fill([1:length(yMean) length(yMean):-1:1],conf, [.5 0.5 0.5],'EdgeColor',[.5 0.5 0.5],'FaceAlpha',.5,'HandleVisibility','off');
        %[h,p,ci,stats] = ttest2(mean(stimresp_c(:,:,2,stimg,stimc),2),mean(avgresp(:,:,2,stimg,stimc),2));
%         for t = 1:length(yMean)     % CoCoSys Lab
%             [~, pval(t)] = ttest(squeeze(cat_avgresp(t,:,stimg,stimc)), squeeze(fix_avgresp(t,:,stimg,stimc)));
%         end
%         % convert to logical
%         signific = nan(1, length(yMean)); 
%         signific(pval < 0.05) = 1;
%         plot(signific * -0.03, '.k');
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
    
    for stimc = 1:length(stimgroups{stimg})
        
        subplot(2,length(stimgroups{stimg}), length(stimgroups{stimg})+stimc); hold on;
        ax = axis;
        axis([ax(1:2) mn mx]);
        %exportgraphics(gca,[sprintf('FvsC/%s%s/',stimgrnames{stimg},stimleg{stimg}(stimc)) 'SF.eps']);

    end
    %figurewrite(sprintf('%s',stimgrnames{stimg}),[],[],plotname);
end

