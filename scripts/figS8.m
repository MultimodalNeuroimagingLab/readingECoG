%% Temporal characteristics of LFB (Fig. S8) 
% 
% time-to-onset, time-to-peak and response duration during fixaiton vs categorizaiton
%
% ZQ 2024 Mayo Clinic
%
%%
close all

winlen_tc = find(epochtime_bb >= -0.5 & epochtime_bb <= 2.5);   % -0.5 s to 2.5 s
winlen_pt = find(epochtime_bb >= 0.25 & epochtime_bb <= 1.25);      % winlen for peak time

% Which ROI elec to select (elec_EVC, elec_VTC, elec_VTC_f, elec_VTC_w, elec_VTC_ns)
% ROI Names:               (     EVC,      VTC,        FFA,       vWFA,       nsVTC)

elec_consider = elec_VTC_w; 
roiName = 'word';

% these are the stimuli we care about
stim = 4:19;

% Combining channels and trials within ROI and grouping stim in categories

%  time x number of visual-responsive electrodes in ROI x tasks x trials x categories x stim
avgresp = nan(length(winlen_pt), length(elec_consider{1})+length(elec_consider{2}), numtasks, numreps, 16);

counter = 1;
for zzz = 1 : 2
    for ccc = elec_consider{zzz}
        avgresp(:, counter, :, :, :) = ...
            squeeze(lfbdata_pc(winlen_pt, zzz, ccc, :, :, stim));
        counter = counter+1;
    end
end
 
assert(counter-1 == length(elec_consider{1})+length(elec_consider{2}));

% el X task X stim X [duration, on_time, pk_time]
met_data = nan(length(elec_consider{1})+length(elec_consider{2}), 2, 16, 3);

% loop across all responsive channels within this ROI
for temp_el = 1 : length(elec_consider{1})+length(elec_consider{2})

    % loop across tasks
    for temp_tsk = 1:2
    
        % loop across all stim
        for stimc = 1 : length(stim)

            temp_data = squeeze(avgresp(:, temp_el, temp_tsk, :, stimc));
            
            temp_resp = temp_data;

            % removing all nan trials
            temp_resp( :, find( all( isnan(temp_data), 1))) = [];
    
            % timepoints for response significantly different from baseline
           % one-sided cluster-based permutation testing (<0, negative one-sided)
            clusters = ieeg_clusterTest1( temp_resp, 0.05, 0.05, 10000, -1); 
    
            temp_dur = 0;
            mn = 0;
            temp_m = squeeze(geomean( temp_resp+1, 2, 'omitnan') -1);

            % loop across significant clusters
            % duration is total length of all significant clusters
            % T2Pk is time of min response within all significant clusters
            for tempkk = 1: size(clusters, 1)
                temp_dur = temp_dur + clusters(tempkk,2) - clusters(tempkk,1);
                if min(temp_m( clusters(tempkk,1): clusters(tempkk,2))) < mn
                    [mn, mn_idx] = min(temp_m( clusters(tempkk,1): clusters(tempkk,2)));
                    temp_pk = clusters(tempkk,1)+mn_idx-1;
                end
            end

            if ~isempty(clusters)

                % T2On is the onset of first significant cluster
                % store all in met_data for all electrodes within ROI               
                met_data(temp_el, temp_tsk, stimc, :) = [temp_dur epochtime_bb( winlen_pt( [min(clusters(:,1)) temp_pk]))]';

            end

        end
       
    end

end

clear temp*

if matches( roiName, 'word'), met_data_w = met_data;
elseif matches( roiName, 'face'), met_data_f = met_data;
elseif matches( roiName, 'EVC'), met_data_e = met_data;
end

clear met_data


%% Combining faces and words for EVC

stimgroups  = {[ 6 10  7 11   8 12   9 13] [14 17 15 18  16 19  4 5  ]};
stimleg     = {["0%", "25%", "50%", "75%"] ["3%", "5%",  "8%", "100%"]};
stimgrnames = {'PhaseCoherence'                      'Contrast'       };

% stim x task x [duration, on_time, pk_time]
met_mean = nan(8, 2, 3);
met_sem = nan(8, 2, 3);
counter = 1;

for stimg = 1 : length(stimgroups)

    for stimc = 1 : 2: length(stimgroups{stimg})

        temp_data = squeeze( mean( met_data_e(:, :, stimgroups{stimg}(stimc:stimc+1)-3, :), 3, 'omitnan'));

        nTrials = size( temp_data, 1);
        ts = tinv([0.16, 0.84], nTrials - 1);

        temp_mean = squeeze( median( temp_data, 1, 'omitnan'));
        temp_std  = squeeze( std( temp_data, 0, 1, "omitnan") / sqrt(nTrials));

        met_mean( counter, :, :) = (temp_mean);
        met_sem ( counter, :, :)  = (ts(1)*temp_std);

        counter = counter+1;

    end

end

assert(counter-1 == 8);

clear temp*

%% Combining faces (for face-selective VTC) and words (for word-selective VTC)
% runthi ssection after running the first section for both face- and word-selective VTC

stimgroups  = {[ 6      7      8      9  ] [10     11     12     13  ] [14    15     16      4   ] [17    18     19      5   ]};
stimleg     = {["0%", "25%", "50%", "75%"] ["0%", "25%", "50%", "75%"] ["3%", "5%",  "8%", "100%"] ["4%", "6%", "10%", "100%"]};
stimgrnames = {'WordPhaseCoherence'        'FacePhaseCoherence'               'WordContrast'              'FaceContrast'};

% stim x task x [duration, on_time, pk_time]
met_mean = nan(8, 2, 3);
met_sem = nan(8, 2, 3);
counter = 1;

for stimg = 1 : 2: length(stimgroups)

    for stimc = 1 : length(stimgroups{stimg})

        temp_data_w = ( squeeze( met_data_w(:, :, stimgroups{stimg}(stimc)-3, :)));
        temp_data_f = ( squeeze( met_data_f(:, :, stimgroups{stimg+1}(stimc)-3, :)));
        temp_data = cat( 1, temp_data_w, temp_data_f);

        nTrials = size( temp_data, 1);
        ts = tinv([0.16, 0.84], nTrials - 1);

        temp_mean = squeeze( median( temp_data, 1, 'omitnan'));
        temp_std  = squeeze( std( temp_data, 0, 1, "omitnan") / sqrt(nTrials));

        met_mean( counter, :, :) = (temp_mean);
        met_sem ( counter, :, :)  = (ts(1)*temp_std);

        counter = counter+1;

    end

end

assert(counter-1 == 8);

clear temp*



%% Plot the temporal characteristics

% temp_dir = fullfile( dataPath.output, 'eps/timings', roiName);
% if not(isfolder(temp_dir))
%     mkdir(temp_dir)
% end

% contrast
stim = 5:8;

% Phase coherence
%stim = [1:4 8];


for ii = 1:3

    % converison factor for duration, on_time, pk_time in ms

    if ii == 1
        scl = 5;
    else
        scl=1000;
    end

    figure

    tempy = squeeze( met_mean( stim, :, ii)) * scl;
    
    b = bar(tempy);
    
    b(1).EdgeColor = "none";
    b(1).FaceColor = [1 0 0];
    b(2).EdgeColor = "none";
    b(2).FaceColor = [0 1 0];


    hold on

    tempx = nan(size(tempy,2), size(tempy,1));
    for jj = 1:size(tempy, 2)
        tempx(jj,:) = b(jj).XEndPoints;
    end

    % Plot the errorbars
    errorbar(tempx', tempy, squeeze( met_sem( stim, :, ii)) * scl, 'k', 'LineStyle', 'none', 'LineWidth', 2, CapSize=0);

    legend('Fix', 'Cat', '', '')

    ylim([0 800])
    ylabel('Time (ms)')

    switch ii

        case 1
            title('Duration')
            % exportgraphics(gca,[temp_dir '/duration.eps'], 'BackgroundColor','none','ContentType','vector', 'Resolution',300);
        case 2
            title('Time-to-onset')
            % exportgraphics(gca,[temp_dir '/on_time.eps'], 'BackgroundColor','none','ContentType','vector', 'Resolution',300);
        case 3
            title('Time-to-peak')
            % exportgraphics(gca,[temp_dir '/pk_time.eps'], 'BackgroundColor','none','ContentType','vector', 'Resolution',300);
        otherwise
            assert(1==0, 'ii should be {1,2,3}')
    end

end


%% Stats
% linear mixed-effects model that predicted time-to-onset, time-to-peak, and duration, independently
% fixed effects for stimulus-contrast X task X ROI (EVC/VTC) 
% and an intercept-only random effect for subjects

%% EVC

stimgroups  = {[ 6 10  7 11   8 12   9 13] [14 17 15 18  16 19  4 5  ]};
stimleg     = {["0%", "25%", "50%", "75%"] ["3%", "5%",  "8%", "100%"]};
stimgrnames = {'PhaseCoherence'                      'Contrast'       };

data_EVC = nan(4, 12, 2, 3);

stimg = 2;

    for stimc = 1 : 2: length(stimgroups{stimg})

        temp_data = squeeze( mean( met_data_e(:, :, stimgroups{stimg}(stimc:stimc+1)-3, :), 3, 'omitnan'));

        data_EVC ((stimc+1)/2, :, :, :) = temp_data;

    end

%% VTC

stimgroups  = {[ 6      7      8      9  ] [10     11     12     13  ] [14    15     16      4   ] [17    18     19      5   ]};
stimleg     = {["0%", "25%", "50%", "75%"] ["0%", "25%", "50%", "75%"] ["3%", "5%",  "8%", "100%"] ["4%", "6%", "10%", "100%"]};
stimgrnames = {'WordPhaseCoherence'        'FacePhaseCoherence'               'WordContrast'              'FaceContrast'};

data_VTC = nan(4, 14, 2, 3);

stimg = 3;

for stimc = 1 : length(stimgroups{stimg})

        temp_data_w = ( squeeze( met_data_w(:, :, stimgroups{stimg}(stimc)-3, :)));
        temp_data_f = ( squeeze( met_data_f(:, :, stimgroups{stimg+1}(stimc)-3, :)));
        temp_data = cat( 1, temp_data_w, temp_data_f);

        data_VTC (stimc, :, :, :) = temp_data;

end

%% Setting up the data matrix

data_timings = cat(2, data_EVC, data_VTC);

for ii = 1:3
    % duration, on_time, pk_time in ms

    if ii == 1
        scl = 5;
    else
        scl=1000;
    end

    data_timings(:,:,:,ii) = scl*data_timings(:,:,:,ii);


    clear ii scl temp*

    % ME stats
    %  4 (contrast levels) x 20 (electrodes) x 2 (fix/cat) x 3 (dur, on, peak)
    yy = squeeze(data_timings(:,:,:,ii));
    [nCon, nElec, nTasks] = size(yy);

    % electrodes wrt subjects (1/2)
    elecSub = [0 0 0 0 1 1 0 0 1 1 1 1 1 0 0 0 0 0 0 0];

    % electrodes wrt ROIs (EVC/VTC)
    elecROI = [zeros(1,6), ones(1,14)];

    % contrast levels (0 1 2 3)
    cond_con   = repmat( reshape(0:(nCon-1), nCon, 1, 1),   1, nElec, nTasks) - 1.5;

    % task (Fix/Cat)
    cond_tasks = repmat( reshape(0:(nTasks-1), 1, 1, nTasks), nCon, nElec, 1);

    % ROI (EVC/VTC)
    cond_ROI   = repmat( reshape(elecROI, 1, nElec, 1),      nCon, 1, nTasks);

    % Subject (1/2)
    cond_sub   = repmat( reshape(elecSub, 1, nElec, 1),      nCon, 1, nTasks);

    % electodes as conditional
    cond_elec  = repmat(  reshape(0:(nElec-1), 1, nElec, 1), nCon, 1, nTasks);

    tbl = table( ...
        yy(:), ...
        cond_con(:), ...
        categorical(cond_tasks(:), 0:1, {'Fix', 'Cat'}), ...
        categorical(cond_ROI(:),   0:1, {'EVC', 'VTC'}), ...
        categorical(cond_sub(:),   0:1, {'Sub-01','Sub-02'}), ...
        categorical(cond_elec(:)), ...
        'VariableNames', ...
        {'Response','Contrast','Tasks','ROI','Subject','Electrode'});

    lme = fitlme(tbl,   ...
        'Response ~ Contrast*Tasks*ROI + (1|Subject) + (1|Electrode)', ...
        'FitMethod','REML', ...
        'DummyVarCoding','effects')

end
