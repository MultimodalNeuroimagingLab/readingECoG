%% PSC-based stimulus and task effects in ECoG HFB/LFB (Fig. S6c,d, S9) 
% 
% time-average HFB/LFB power changes computed from power spectral components during fixaiton vs categorizaiton
% this can be run as an independent script
%
% ZQ 2023 Mayo Clinic
%
%%
% We compute psd separately in 2 temporal windows (0:1:2 s)
% This is an independent script so setup and load data

clc; clear; close all;

dataPath = setLocalDataPath(1);
load(fullfile( dataPath.input, 'derivatives/psd_timechunks1s.mat'))

%% HFB & LFB analysis at different time windows using spectral-PCA

% gcc = {[72, 74, 75, 76], [106, 108]};   % primary visual
% roiName = 'EVC'; 

% gcc =  {[10 12 41 45 46 77 78],[]};     % face-VTC
% roiName = 'face';

gcc =  {[11 44],[20 22 84 86 110]};   % word-VTC
roiName = 'word';

avgresp_bb = nan(length(gcc{1})+length(gcc{2}), numreps, numtasks, 16);    % 16 = # of stim 4-19
avgresp_nb = nan(length(gcc{1})+length(gcc{2}), numreps, numtasks, 16);

% mark bad trials
bd_tr = [ 30,  54, 136, ... 
         222, 246, 328];

counter = 1;

% loop across subjects
for zzz = 1:length(gcc)

    % loop across channels in this ROI
    for ccc = gcc{zzz} 

        fprintf('\nSub:%d, Ch:%03d',zzz,ccc)

        test_psd = [reshape( psd_on( zzz, ccc, :, :, 4:19, 6:200, 1), [], 195); 
                    reshape( psd_off(zzz, ccc, :, :, 4:19, 6:200), [], 195)];

        test_psd(bd_tr, :) = nan;

        avg_psd = mean(test_psd, 1, 'omitnan');

        % log-normalize
        psd_norm = log10( test_psd) - log10( avg_psd);

        % PCA
        [coeff,score,latent,tsquared,explained] = pca(psd_norm, 'Centered', false);

        % Remove bands with notch
        coeff(50:70, :) = 0;
        coeff(112:124, :) = 0;
        coeff(172:192, :) = 0;

        % First PC - HFB; Second PC- LFB
        psc1 = coeff(:,1);
        psc2 = coeff(:,2);

        % Plot

        % figure,
        % hold on
        % for jj = 1:2
        %     plot(f_psd(6:200), coeff(:,jj))
        % end
        % yline(0, 'k--')
        % 
        % xlabel('Frequency(Hz)')
        % ylabel('log-Power')
        % title(sprintf('Sub:%d, Ch:%03d',zzz,ccc))
        % 
        % legend(sprintf('First PC (%0.2f)', explained(1)), sprintf('Second PC (%0.2f)', explained(2)))

        wt_bb = psd_norm * psc1;
        wt_nb = psd_norm * psc2;

        wt_bb(bd_tr) = NaN;
        wt_nb(bd_tr) = NaN;

        % only stim trials
        xyz = reshape(wt_bb(1:192), 2, 6, 16, []);

        avgresp_bb(counter, :, :, :) = permute(xyz, [2 1 3]);

        xyz = reshape(wt_nb(1:192), 2, 6, 16, []);

        avgresp_nb(counter, :, :, :) = permute(xyz, [2 1 3]);

        counter = counter+1;

    end

end

assert(counter-1 == length(gcc{1})+length(gcc{2}));


%% HFB analysis

%  Stimulus Categories

% For contrast
% stimgroups  = {[14    15     16      4   ] [17    18     19      5   ]};
% stimleg     = {["3%", "5%",  "8%", "100%"] ["4%", "6%", "10%", "100%"]};
% stimgrnames = {  'WordContrast'              'FaceContrast'};


% % For PC
stimgroups  = {[ 6      7      8      9       4  ] [10     11     12     13       5  ] };
stimleg     = {["0%", "25%", "50%", "75%", "100%"] ["0%", "25%", "50%", "75%", "100%"] };
stimgrnames = {'WordPhaseCoherence'               'FacePhaseCoherence'             };


% temp_dir = fullfile( dataPath.output, 'eps/avg/psc_bb', roiName, winlab);
% if not(isfolder(temp_dir))
%     mkdir(temp_dir)
% end


avgresp = normalize( reshape (avgresp_bb, [], 1), 'range');

avgresp = reshape (avgresp, [], 2, 16);

for stimg = 1:length(stimgroups)

    temp_yMeanc =[];
    temp_yMeanf =[];
    temp_ySEMc =[];
    temp_ySEMf =[];
    temp_posx  =[];
    posx = 0;
    figure

    stim_ix = stimgroups{stimg};
    
    for stimc = 1 : length(stim_ix)

        % Fixation

        tempData = avgresp(:, 1, stim_ix(stimc)-3);
        
        nTrials = size( tempData, 1) - sum( isnan( tempData));
        
        % Median, SEM and 68% CI
        yMean = median( tempData, 1, 'omitnan') ;
        ySEM  =  std( tempData, 0, 1 , 'omitnan') / sqrt(nTrials);
        ts = tinv([0.16, 0.84], nTrials - 1);

        temp_yMeanf =[temp_yMeanf yMean];
        temp_ySEMf  =[temp_ySEMf  (ts(1)*ySEM)'];
      
        % Categorization

        tempData = avgresp(:, 2, stim_ix(stimc)-3);
        
        nTrials = size( tempData, 1) - sum( isnan( tempData));
        
        % Median, SEM and 68% CI
        yMean = median( tempData, 1, 'omitnan') ;
        ySEM  =  std( tempData, 0, 1 , 'omitnan') / sqrt(nTrials);
        ts = tinv([0.16, 0.84], nTrials - 1);

        temp_yMeanc =[temp_yMeanc yMean];
        temp_ySEMc  =[temp_ySEMc  (ts(1)*ySEM)'];
        
        posx = posx + 1;

        temp_posx = [temp_posx posx];
        
    end

    hold on;
    errorbar(temp_posx, temp_yMeanf, temp_ySEMf, Color = [1 0 0], LineStyle = '-', CapSize=0);
    errorbar(temp_posx, temp_yMeanc, temp_ySEMc, Color = [0 1 0], LineStyle = '-', CapSize=0);
    

    xlim([0.5 stimc+0.5])
    ylim([0.2 0.7])   % EVC Con
    %ylim([0.15 .75])
    set(gca,'XTick',1:stimc);
    xticklabels(stimleg{stimg})
    ylabel('PSC HFB power (x-fold increase)')
    xlabel(stimgrnames{stimg})

    % exportgraphics(gca,[temp_dir sprintf('/avg%s.eps', stimgrnames{stimg})], 'BackgroundColor','none','ContentType','vector', 'Resolution',300);

    % Scaling Factor
    
    figure, 
    plot((temp_yMeanc)./(temp_yMeanf), 'c-');
    xlim([0.5 stimc+0.5])
    ylim([0.9 1.6])
    yline(1, 'k--')
    set(gca,'XTick',1:stimc);
    xticklabels(stimleg{stimg})
    ylabel('Scaling Factor')
    xlabel(stimgrnames{stimg})

    % exportgraphics(gca,[temp_dir sprintf('/SF%s.eps', stimgrnames{stimg})], 'BackgroundColor','none','ContentType','vector', 'Resolution',300);

end



clear temp*

%% LFB analysis
% similar to HFB analysis

%  Stimulus Categories

% For contrast
stimgroups  = {[14    15     16      4   ] [17    18     19      5   ]};
stimleg     = {["3%", "5%",  "8%", "100%"] ["4%", "6%", "10%", "100%"]};
stimgrnames = {  'WordContrast'              'FaceContrast'};


% For PC
% stimgroups  = {[ 6      7      8      9       4  ] [10     11     12     13       5  ] };
% stimleg     = {["0%", "25%", "50%", "75%", "100%"] ["0%", "25%", "50%", "75%", "100%"] };
% stimgrnames = {'WordPhaseCoherence'               'FacePhaseCoherence'             };


% temp_dir = fullfile( dataPath.output, 'eps/avg/psc_nb', roiName, winlab);
% if not(isfolder(temp_dir))
%     mkdir(temp_dir)
% end

avgresp = normalize( reshape (avgresp_nb, [], 1), 'range', [-1 0]);
avgresp = reshape (avgresp, [], 2, 16);

for stimg = 1:length(stimgroups)

    temp_yMeanc =[];
    temp_yMeanf =[];
    temp_ySEMc =[];
    temp_ySEMf =[];
    temp_posx  =[];
    posx = 0;
    figure

    stim_ix = stimgroups{stimg};
    
    for stimc = 1 : length(stim_ix)

        % Fixation

        tempData = avgresp(:, 1, stim_ix(stimc)-3);
        
        nTrials = size( tempData, 1) - sum( isnan( tempData));
        
        % Median, SEM and 68% CI
        yMean = median( tempData, 1, 'omitnan') ;
        ySEM  =  std( tempData, 0, 1 , 'omitnan') / sqrt(nTrials);
        ts = tinv([0.16, 0.84], nTrials - 1);

        temp_yMeanf =[temp_yMeanf yMean];
        temp_ySEMf  =[temp_ySEMf  (ts(1)*ySEM)'];
      
        % Categorization

        tempData = avgresp(:, 2, stim_ix(stimc)-3);
        
        nTrials = size( tempData, 1) - sum( isnan( tempData));
        
        % Median, SEM and 68% CI
        yMean = median( tempData, 1, 'omitnan') ;
        ySEM  =  std( tempData, 0, 1 , 'omitnan') / sqrt(nTrials);
        ts = tinv([0.16, 0.84], nTrials - 1);

        temp_yMeanc =[temp_yMeanc yMean];
        temp_ySEMc  =[temp_ySEMc  (ts(1)*ySEM)'];
        
        posx = posx + 1;

        temp_posx = [temp_posx posx];
        
    end

    hold on;
    errorbar(temp_posx, temp_yMeanf, temp_ySEMf, Color = [1 0 0], LineStyle = '-', CapSize=0);
    errorbar(temp_posx, temp_yMeanc, temp_ySEMc, Color = [0 1 0], LineStyle = '-', CapSize=0);
    
    set(gca,'Ydir','reverse')

    xlim([0.5 stimc+0.5])
    ylim([-.75 -.35])
    yline(0, 'k--')
    set(gca,'XTick',1:stimc);
    xticklabels(stimleg{stimg})
    ylabel('PSC LFB power (x-fold decrease)')
    xlabel(stimgrnames{stimg})

    % exportgraphics(gca,[temp_dir sprintf('/avg%s.eps', stimgrnames{stimg})], 'BackgroundColor','none','ContentType','vector', 'Resolution',300);

end



clear temp*

