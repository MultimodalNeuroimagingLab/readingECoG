%% Stimulus and task effects in ECoG LFB (Fig. S7) 
% 
% time resolved LFB power changes during fixaiton vs categorizaiton
%
% ZQ 2023 Mayo Clinic
%
%% average across all visual-responsive channels

close all

% Which ROI elec to select (elec_EVC, elec_VTC, elec_VTC_f, elec_VTC_w, elec_VTC_ns, elec_ROI)
% ROI Names:               (     EVC,      VTC,        FFA,       vWFA,       nsVTC, allResp)

elec_consider = elec_ROI; 
roiName = 'allResp';


%  Stimulus Categories


stimgroups  = {[ 6      7      8      9       4  ] [10     11     12     13       5  ] [14    15     16      4   ] [17    18     19      5   ]};
stimleg     = {["0%", "25%", "50%", "75%", "100%"] ["0%", "25%", "50%", "75%", "100%"] ["3%", "5%",  "8%", "100%"] ["4%", "6%", "10%", "100%"]};
stimgrnames = {'WordPhaseCoherence'               'FacePhaseCoherence'               'WordContrast'              'FaceContrast'};


winlen_tc = find(epochtime_bb >= -0.5 & epochtime_bb <= 2.5);   % -0.5 s to 2.5 s
winlen_sf = find(epochtime_bb >= 0.25 & epochtime_bb <= 1.25);        %    0 s to 1 s

% Combining channels and trials within ROI and grouping stim in categories

%  time x number of visual-responsive electrodes in ROI x tasks x trials x categories x stim
avgresp = nan(length(winlen_tc), length(elec_consider{1})+length(elec_consider{2}), numtasks, numreps, length(stimgroups), 5);
peaktime = nan(numtasks+2, length(stimgroups), 5);

counter = 1;
for zzz = 1 : 2
    for ccc = elec_consider{zzz}
        for stimg = 1 : length(stimgroups)
            for stimc = 1 : length(stimgroups{stimg})
                avgresp(:, counter, :, :, stimg, stimc) = ...
                    squeeze(lfbdata_pc(winlen_tc, zzz, ccc, :, :, stimgroups{stimg}(stimc)));
            end
        end
        counter = counter+1;
    end
end
 
assert(counter-1 == length(elec_consider{1})+length(elec_consider{2}));

%  merge all trials, so reshape as time x [responsive_channels x trials] x task x categories x stim
avgresp = permute(avgresp, [1 2 4 3 5 6]);
avgresp = reshape(avgresp, [length(winlen_tc), (length(elec_consider{1})+length(elec_consider{2})) * numreps, numtasks, length(stimgroups), 5]);

% Plot

% smoothing window (1 = no smoothing)
smooth_len = 1;

% loop across stimulus groups
for stimg = 3:4%1 : length(stimgroups)

    % temp_dir = fullfile( dataPath.output, 'eps/nbTC', roiName, stimgrnames{stimg});
    % if not(isfolder(temp_dir))
    %     mkdir(temp_dir)
    % end

    figure('Position', [0 0 250*length(stimgroups{stimg}) 500], 'Visible', 'on');
    
    mx = -Inf; mn = Inf;
    y_fix = [];
    y_cat = [];

    % for stimuli within the selected stimulus group
    for stimc = 1 : length(stimgroups{stimg})
        
        subplot(2,length(stimgroups{stimg}),stimc); hold on;

        % Fixation
        
        data_smooth = log(movmean( squeeze( avgresp( :, :, 1, stimg, stimc)), smooth_len, 1, "omitnan") + 1);

        nTrials = size( data_smooth, 2) - sum( all( isnan(data_smooth), 1));
        
        % Mean, SEM and 68% CI
        yMean = mean( data_smooth, 2, "omitnan");
        y_fix = [y_fix yMean(winlen_pt)];
        ySEM  = std( data_smooth, 0, 2, "omitnan") / sqrt(nTrials);
        ts = tinv([0.16, 0.84], nTrials - 1);
        conf = exp([(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'])-1;

        % plot in linear scale
        plot(exp(yMean)-1, 'r-', 'LineWidth', 1);
        fill([1:length(winlen_tc) length(winlen_tc):-1:1],conf, [.5 0 0],'EdgeColor','none','FaceAlpha',.5,'HandleVisibility','off');
        
        % Categorization

        data_smooth = log(movmean( squeeze( avgresp( :, :, 2, stimg, stimc)), smooth_len, 1, "omitnan") + 1);

        nTrials = size( data_smooth, 2) - sum( all( isnan(data_smooth), 1));
        
        % Mean, SEM and 68% CI
        yMean = mean( data_smooth, 2, "omitnan");
        y_cat = [y_cat yMean(winlen_pt)];
        ySEM  = std( data_smooth, 0, 2, "omitnan") / sqrt(nTrials);
        ts = tinv([0.16, 0.84], nTrials - 1);
        conf = exp([(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'])-1;

        % plot in linear scale
        plot(exp(yMean)-1, 'g-', 'LineWidth', 1);
        fill([1:length(winlen_tc) length(winlen_tc):-1:1],conf, [0 0.5 0],'EdgeColor','none','FaceAlpha',.5,'HandleVisibility','off');
        
        yline(0, '--k')
        xline(find(epochtime_bb == 0),'--k');
        ax = axis;
        mx = max(mx,max(ax(3:4)));
        mn = min(mn,min(ax(3:4)));
        xlim([1 length(winlen_tc)]);
        ylabel('HFB Power (x-fold increase)');
        xlabel('Time (s)');
        set(gca,'XTick',1:100:length(winlen_tc)); set(gca,'XTickLabel',epochtime_bb(winlen_tc(1)):0.5:epochtime_bb(winlen_tc(end)));
        title(stimgrnames{stimg}+" "+stimleg{stimg}(stimc));
        legend('F', 'C', '', '');
        
    end
    
     % adjust y-range for all plots
    for stimc = 1:length(stimgroups{stimg})
        
        subplot(2,length(stimgroups{stimg}),stimc); hold on;
        ax = axis;
        axis([ax(1:2) -0.6 0.6]);
        % exportgraphics(gca,[temp_dir sprintf('/STC%s%s.eps',stimgrnames{stimg},stimleg{stimg}(stimc))], 'BackgroundColor','none','ContentType','vector','Resolution',300);
    
    end
        
    mx = -Inf; mn = Inf;
    
    % for stimuli within the selected stimulus group
    for stimc = 1 : length(stimgroups{stimg})
        
        subplot(2,length(stimgroups{stimg}), length(stimgroups{stimg})+stimc); hold on;

        % Scaling Factor

        scaling_factor = [];

        
        tempc = squeeze( avgresp(winlen_sf, :, 2, stimg, stimc));
        tempf = squeeze( avgresp(winlen_sf, :, 1, stimg, stimc)); 

        yMeanc = movmean( geomean( tempc+1, 2, 'omitnan')-1, 40, 1, "omitnan");
        yMeanf = movmean( geomean( tempf+1, 2, 'omitnan')-1, 40, 1, "omitnan");
        % 
        % % for cases where fix resp < 0, in FFA & nsVTC for 3% WC
        % yMeanf( find( yMeanf < 0.01)) = 0.01;
        % 
        % if(sum(yMeanf < 0.01)>0)
        %     fprintf('\n%d,%d', stimg, stimc)
        % end

        yMeansf = ones(size(yMeanc));

        tempsc = movmean( zerodiv( yMeanc, yMeanf, 1), 1, 1, "omitnan");

        % Remove trials that are all Nan
        tempf( :, find( all( isnan(tempf), 1))) = [];
        tempc( :, find( all( isnan(tempc), 1))) = [];

        clusters = ieeg_clusterTest(tempf, tempc); 

       % cluster-based permutation testing to find instance where categorizaiton responses significantly differ from fixaiton
        for tempkk = 1 : size(clusters, 1)
            yMeansf(clusters(tempkk,1):clusters(tempkk,2)) = tempsc(clusters(tempkk,1):clusters(tempkk,2));  
        end


        % if ~isempty(clusters)
        %     t_sf = [t_sf; [stimgroups{stimg}(stimc) epochtime_bb( winlen_sf( clusters(1,1)))]];
        % end

        plot( yMeansf, 'c-', 'LineWidth', 1);

        yline(1, '--k')
        ax = axis;
        axis([100 300 ax(3:4)]);
        mx = max(mx,max(ax(3:4)));
        mn = min(mn,min(ax(3:4)));
        xlim([1 length(winlen_sf)]);
        ylabel('Scaling Factor');
        xlabel('Time (s)');
        set(gca,'XTick',1:100:length(winlen_sf)); set(gca,'XTickLabel',epochtime_bb(winlen_sf(1)):0.5:epochtime_bb(winlen_sf(end)));
        
    end
    
    % adjust y-range for all plots
    for stimc = 1:length(stimgroups{stimg})
        
        subplot(2,length(stimgroups{stimg}), length(stimgroups{stimg})+stimc); hold on;
        ax = axis;
        axis([ax(1:2) 0.5 3]);
        %exportgraphics(gca,[temp_dir sprintf('/SF%s%s.eps',stimgrnames{stimg},stimleg{stimg}(stimc))], 'BackgroundColor','none','ContentType','vector', 'Resolution',300);

    end

 
end

clear temp*



