%% Stimulus and task effects in ECoG HFB/LFB and in fMRI BOLD (Fig. 5, 6, S6 a,b, S9) 
% 
% time-average HFB/LFB power changes during fixaiton vs categorizaiton
%
% ZQ 2023 Mayo Clinic
%
%% fMRI BOLD Results - Avg. responses per ROI

% this section reproduces results from
% Kendrick N Kay & Jason D Yeatman (2017) Bottom-up and top-down computations in word- and face-selective cortex eLife 6:e22341.

load( fullfile( dataPath.input, 'derivatives', 'BOLD.mat'))

%ecog_orderWRTfMRI = [1 17 18 19 5 14 15 16 4  20 21 22 10 11 12 13 6 7 8 9 23 2 24];

close all

% Which ROI to select (1:4,   6,    5,   8)
% ROI Names:          (EVC, FFA, vWFA, IPS)

roi_consider = 5;
roiName = 'vWFA';

%  Stimulus Categories

% For Contrast
stimgroups  = {[6      7     8      9    ] [2     3     4      5     ]};
stimleg     = {["3%", "5%",  "8%", "100%"] ["4%", "6%", "10%", "100%"]};
stimgrnames = {'WordContrast'              'FaceContrast'};

% For Phase Coherence
% stimgroups  = {[ 17    18      19    20       9  ] [13     14     15     16       5  ]};
% stimleg     = {["0%", "25%", "50%", "75%", "100%"] ["0%", "25%", "50%", "75%", "100%"]};
% stimgrnames = {'WordPhaseCoherence'               'FacePhaseCoherence'};


temp_dir = fullfile( dataPath.output, 'eps/avg/fmri', roiName);
if not(isfolder(temp_dir))
    mkdir(temp_dir)
end
 
% loop across stimulus groups
for stimg = 1:length(stimgroups)

    figure

    temp_yMeanc = [];
    temp_yMeanf = [];
    temp_ySEMc  = [];
    temp_ySEMf  = [];
    temp_posx   = [];
    posx = 0;

    % for stimuli within the selected stimulus group
    for stimc = 1 : length(stimgroups{stimg})

        % Fixation - Mean & SEM
        temp_yMeanf = [temp_yMeanf mean( squeeze(   groupbeta( roi_consider, stimgroups{stimg}(stimc), 1)), 1)];
        temp_ySEMf  = [temp_ySEMf  mean( squeeze( groupbetase( roi_consider, stimgroups{stimg}(stimc), 1)), 1)];

        % Categorization - Mean & SEM
        temp_yMeanc = [temp_yMeanc mean( squeeze(   groupbeta( roi_consider, stimgroups{stimg}(stimc), 2)), 1)];
        temp_ySEMc  = [temp_ySEMc  mean( squeeze( groupbetase( roi_consider, stimgroups{stimg}(stimc), 2)), 1)];
        
        % x-axis locations
        posx = posx + 1;
        temp_posx  = [temp_posx posx];
        
    end


    % plot
    hold on
    errorbar(temp_posx, temp_yMeanf, temp_ySEMf, 'r-', CapSize=0);
    errorbar(temp_posx, temp_yMeanc, temp_ySEMc, 'g-', CapSize=0);
    xlim([0.5 stimc+0.5])
    ylim([0 1.8])
    set(gca,'XTick',1:stimc);
    xticklabels(stimleg{stimg})
    ylabel(['fMRI BOLD (% change)'])
    xlabel(stimgrnames{stimg})

    % exportgraphics(gca,[temp_dir sprintf('/avg%s.eps',stimgrnames{stimg})], 'BackgroundColor','none','ContentType','vector', 'Resolution',300);

    % Scaling Factor
    
    figure, 
    plot((temp_yMeanc)./(temp_yMeanf), 'b-');
    xlim([0.5 stimc+0.5])
    ylim([0.5 3.5])
    set(gca,'XTick',1:stimc);
    xticklabels(stimleg{stimg})
    ylabel('Scaling Factor')
    xlabel(stimgrnames{stimg})
    yline(1, 'k--')

    % exportgraphics(gca,[temp_dir sprintf('/SF%s.eps',stimgrnames{stimg})], 'BackgroundColor','none','ContentType','vector', 'Resolution',300);


end


clear temp*

%% HFB averaged over time window
 
close all

% Which ROI to select (elec_EVC, elec_VTC, elec_VTC_f, elec_VTC_w, elec_VTC_ns)
% ROI Names:          (     EVC,      VTC,       face,       word,       nsVTC)

elec_consider = elec_VTC_w; 
roiName = 'word';


%  Stimulus Categories

% For contrast
% stimgroups  = {[14    15     16      4   ] [17    18     19      5   ]};
% stimleg     = {["3%", "5%",  "8%", "100%"] ["4%", "6%", "10%", "100%"]};
% stimgrnames = {  'WordContrast'              'FaceContrast'};


% For PC
stimgroups  = {[ 6      7      8      9       4  ] [10     11     12     13       5  ] };
stimleg     = {["0%", "25%", "50%", "75%", "100%"] ["0%", "25%", "50%", "75%", "100%"] };
stimgrnames = {'WordPhaseCoherence'               'FacePhaseCoherence'             };


winlen_resp = find(epochtime_bb >= 0 & epochtime_bb <= 1);   

% temp_dir = fullfile( dataPath.output, 'eps/avg/bb', roiName);
% if not(isfolder(temp_dir))
%     mkdir(temp_dir)
% end

% Combining trials and channels

%  time x number of visually responsive channels x tasks x trials x categories x stim

avgresp = nan(length(elec_consider{1})+length(elec_consider{2}), numreps, numtasks, length(stimgroups), 5); 

% Extarcting only channels under consideration over both subjects
counter = 1;
for zzz = 1 : 2
    for ccc = elec_consider{zzz}
        for stimg = 1 : length(stimgroups)
            for stimc = 1 : length(stimgroups{stimg})
                avgresp(counter,  :, :, stimg, stimc) = ...
                    squeeze( geomean( hfbdata_pc( winlen_resp, zzz, ccc, :, :, stimgroups{stimg}(stimc)) + 1, 1, 'omitnan'))';
            end
        end
        counter = counter+1;
    end
end
 
assert(counter-1 == length(elec_consider{1})+length(elec_consider{2}));

% merge all trials within ROI, so reshape as [el_ROI x trials] x task x categories x stim
avgresp = reshape(avgresp, [(length(elec_consider{1})+length(elec_consider{2})) * numreps, numtasks, length(stimgroups), 5]);

% loop across stimulus groups
for stimg = 1:length(stimgroups)

    temp_yMeanc =[];
    temp_yMeanf =[];
    temp_ySEMc =[];
    temp_ySEMf =[];
    temp_posx  =[];
    posx = 0;
    figure
    
    % for stimuli within the selected stimulus group
    for stimc = 1 : length(stimgroups{stimg})

        % Fixation

        % log-transform
        tempData = log( avgresp(:, 1, stimg, stimc));
        
        nTrials = size( tempData, 1) - sum( isnan( tempData));
        
        % Mean, SEM and 68% CI
        yMean = mean( tempData, 1, 'omitnan') ;
        ySEM  =  std( tempData, 0, 1 , 'omitnan') / sqrt(nTrials);
        ts = tinv([0.16, 0.84], nTrials - 1);

        temp_yMeanf =[temp_yMeanf (exp(yMean)-1)];
        temp_ySEMf  =[temp_ySEMf  (exp( yMean+ts(2)*ySEM) - exp(yMean))];
      
        % Categorization

        % log-transform
        tempData = log( avgresp(:, 2, stimg, stimc));
        
        nTrials = size( tempData, 1) - sum( isnan( tempData));
        
        % Mean, SEM and 68% CI
        yMean = mean( tempData, 1, 'omitnan');
        ySEM  =  std( tempData, 0, 1 , 'omitnan') / sqrt(nTrials);
        ts = tinv([0.16, 0.84], nTrials - 1);

        temp_yMeanc =[temp_yMeanc exp(yMean)-1];
        temp_ySEMc  =[temp_ySEMc  (exp( yMean+ts(2)*ySEM) - exp(yMean))];

        
        posx = posx + 1;

        temp_posx = [temp_posx posx];

        
    end
    
    hold on;
    errorbar(temp_posx, temp_yMeanf, temp_ySEMf, Color = [1 0 0], LineStyle = '-', CapSize=0);
    errorbar(temp_posx, temp_yMeanc, temp_ySEMc, Color = [0 1 0], LineStyle = '-', CapSize=0);
    

    xlim([0.5 stimc+0.5])
    ylim([0 3.2])
    set(gca,'YTick',0:1:3);
    %set(gca,'YTick',0:0.4:2.4);
    set(gca,'XTick',1:stimc);
    xticklabels(stimleg{stimg})
    ylabel('ECoG HFB power (x-fold increase)')
    xlabel(stimgrnames{stimg})

    %exportgraphics(gca,[temp_dir sprintf('/avg%s.eps',stimgrnames{stimg})], 'BackgroundColor','none','ContentType','vector', 'Resolution',300);

    % Scaling Factor
    
    figure, 
    plot((temp_yMeanc)./(temp_yMeanf), 'c-');
    xlim([0.5 stimc+0.5])
    ylim([0.5 3])
    yline(1, 'k--')
    set(gca,'XTick',1:stimc);
    xticklabels(stimleg{stimg})
    ylabel('Scaling Factor')
    xlabel(stimgrnames{stimg})

    %exportgraphics(gca,[temp_dir sprintf('/SF%s.eps',stimgrnames{stimg})], 'BackgroundColor','none','ContentType','vector', 'Resolution',300);
    
end

clear temp*


%% LFB averaged over time window
% similar to HFB but with slight offset and reverse y-axis

close all

% Which ROI to select (elec_EVC, elec_VTC, elec_VTC_f, elec_VTC_w, elec_VTC_ns)
% ROI Names:          (     EVC,      VTC,       face,       word,       nsVTC)

elec_consider = elec_VTC_w; 
roiName = 'word';


%  Stimulus Categories

% For contrast
% stimgroups  = {[14    15     16      4   ] [17    18     19      5   ]};
% stimleg     = {["3%", "5%",  "8%", "100%"] ["4%", "6%", "10%", "100%"]};
% stimgrnames = {  'WordContrast'              'FaceContrast'};


% For PC
stimgroups  = {[ 6      7      8      9       4  ] [10     11     12     13       5  ] };
stimleg     = {["0%", "25%", "50%", "75%", "100%"] ["0%", "25%", "50%", "75%", "100%"] };
stimgrnames = {'WordPhaseCoherence'               'FacePhaseCoherence'             };


winlen_resp = find(epochtime_bb >= 0.25 & epochtime_bb <= 1.25);   % 0.25 to 1.25 s, to minimize the effect of evoked responses.

% temp_dir = fullfile( dataPath.output, 'eps/avg/nb', roiName);
% if not(isfolder(temp_dir))
%     mkdir(temp_dir)
% end

% Combining trials and channels

%  time x number of visually responsive channels x tasks x trials x categories x stim

avgresp = nan(length(elec_consider{1})+length(elec_consider{2}), numreps, numtasks, length(stimgroups), 5); 
peaktime = nan(numtasks+2, length(stimgroups), 5);

% Extarcting only channels under consideration over both subjects
counter = 1;
for zzz = 1 : 2
    for ccc = elec_consider{zzz}
        for stimg = 1 : length(stimgroups)
            for stimc = 1 : length(stimgroups{stimg})
                avgresp(counter,  :, :, stimg, stimc) = ...
                    squeeze( geomean( lfbdata_pc( winlen_resp, zzz, ccc, :, :, stimgroups{stimg}(stimc)) + 1, 1, 'omitnan'))';
            end
        end
        counter = counter+1;
    end
end
 
assert(counter-1 == length(elec_consider{1})+length(elec_consider{2}));

avgresp = reshape(avgresp, [(length(elec_consider{1})+length(elec_consider{2})) * numreps, numtasks, length(stimgroups), 5]);


for stimg = 1:length(stimgroups)

    temp_yMeanc =[];
    temp_yMeanf =[];
    temp_ySEMc =[];
    temp_ySEMf =[];
    temp_posx  =[];
    posx = 0;
    figure
    
    for stimc = 1 : length(stimgroups{stimg})

        % Fixation

        tempData = log( avgresp(:, 1, stimg, stimc));
        
        nTrials = size( tempData, 1) - sum( isnan( tempData));

        yMean = mean( tempData, 1, 'omitnan') ;
        ySEM  =  std( tempData, 0, 1 , 'omitnan') / sqrt(nTrials);
        ts = tinv([0.16, 0.84], nTrials - 1);

        temp_yMeanf =[temp_yMeanf exp(yMean)-1];
        temp_ySEMf  =[temp_ySEMf  (exp( yMean+ts(2)*ySEM) - exp(yMean))];

      
        % Categorization

        tempData = log( avgresp(:, 2, stimg, stimc));
        
        nTrials = size( tempData, 1) - sum( isnan( tempData));
        
        yMean = mean( tempData, 1, 'omitnan') ;
        ySEM  =  std( tempData, 0, 1 , 'omitnan') / sqrt(nTrials);
        ts = tinv([0.16, 0.84], nTrials - 1);

        temp_yMeanc =[temp_yMeanc exp(yMean)-1];
        temp_ySEMc  =[temp_ySEMc  (exp( yMean+ts(2)*ySEM) - exp(yMean))];
        
        posx = posx + 1;

        temp_posx = [temp_posx posx];
        
    end

    hold on;
    errorbar(temp_posx, temp_yMeanf, temp_ySEMf, Color = [1 0 0], LineStyle = '-', CapSize=0);
    errorbar(temp_posx, temp_yMeanc, temp_ySEMc, Color = [0 1 0], LineStyle = '-', CapSize=0);
    
    % reverse the y-axis to emphasis decreases (only for visulaization)
    set(gca,'Ydir','reverse')

    xlim([0.5 stimc+0.5])
    ylim([-0.65 0.05])
    yline(0, 'k--')
    set(gca,'XTick',1:stimc);
    xticklabels(stimleg{stimg})
    ylabel('ECoG LFB power (x-fold decrease)')
    xlabel(stimgrnames{stimg})

    %exportgraphics(gca,[temp_dir sprintf('/avg%s.eps',stimgrnames{stimg})], 'BackgroundColor','none','ContentType','vector', 'Resolution',300);

end



clear temp*




%% Fig  6c -- average LFB decreases across all visual-responsive channels and all conditions then do line fit

close all

%  Stimulus Categories

% For contrast
stimgroups  = {[14    15     16      4   ] [17    18     19      5   ]};
stimleg     = {["3%", "5%",  "8%", "100%"] ["4%", "6%", "10%", "100%"]};
stimgrnames = {  'WordContrast'              'FaceContrast'};


winlen_resp = find(epochtime_bb >= 0.25 & epochtime_bb <= 1.25);   % 0.25 to 1.25 s, to minimize the effect of evoked responses.

% temp_dir = fullfile( dataPath.output, 'eps/avg/nb');
% if not(isfolder(temp_dir))
%     mkdir(temp_dir)
% end

% Combining trials and channels

%  time x number of visually responsive channels x tasks x trials x categories x stim

avgresp = nan(length(elec_ROI{1})+length(elec_ROI{2}), numreps, numtasks, length(stimgroups), 4); 

% Extarcting only channels under consideration over both subjects
counter = 1;
for zzz = 1 : 2
    for ccc = elec_ROI{zzz}
        for stimg = 1 : length(stimgroups)
            for stimc = 1 : length(stimgroups{stimg})
                avgresp(counter,  :, :, stimg, stimc) = ...
                    squeeze( geomean( lfbdata_pc( winlen_resp, zzz, ccc, :, :, stimgroups{stimg}(stimc)) + 1, 1, 'omitnan'))';
            end
        end
        counter = counter+1;
    end
end
 
assert(counter-1 == length(elec_ROI{1})+length(elec_ROI{2}));

avgresp = permute(avgresp, [1 2 4 3 5]);

avgresp = reshape(avgresp, [(length(elec_ROI{1})+length(elec_ROI{2})) * numreps * length(stimgroups), numtasks, 4]);

for tsk = 1:2     % Fixation or Categorizaiton

    temp_yMean =[];
    temp_ySEM =[];
    temp_posx  =[];
    posx = 0;
    figure
   
    for stimc = 1 : length(stimgroups{1})

        tempData = log( avgresp(:, tsk, stimc));

        nTrials = size( tempData, 1) - sum( isnan( tempData));

        yMean = mean( tempData, 1, 'omitnan') ;
        ySEM  =  std( tempData, 0, 1 , 'omitnan') / sqrt(nTrials);
        ts = tinv([0.16, 0.84], nTrials - 1);

        temp_yMean =[temp_yMean exp(yMean)-1];
        temp_ySEM  =[temp_ySEM  (exp( yMean+ts(2)*ySEM) - exp(yMean))];

        posx = posx + 1;

        temp_posx = [temp_posx posx];

    end

    % line fit
    y = reshape( squeeze( avgresp(:,tsk,:)), [], 1);
    x = repelem(1:size(avgresp,3), size(avgresp,1))';

    x = x(~isnan(y));
    y = log(y(~isnan(y)));

    [coeffs,S] = polyfit(x, y, 1)
    yfit = polyval(coeffs,x);  

    ste = sqrt(diag(inv(S.R)*inv(S.R')).*S.normr.^2./S.df);

    % Calculate t-statistics
    t_stats = coeffs ./ ste;

    % Calculate p-values using t-distribution
    p_values = 2*(1 - tcdf(abs(t_stats), S.df))

    errorbar(temp_posx, temp_yMean, temp_ySEM, Color = [0 1 0], LineStyle = '-', CapSize=0);
    hold on;
    plot(temp_posx, exp(yfit([find(diff(yfit)); length(yfit)]))-1, 'k-')   

    set(gca,'Ydir','reverse')

    xlim([0.5 stimc+0.5])
    ylim([-0.65 0.05])
    yline(0, 'k--')
    set(gca,'XTick',1:stimc);
    xticklabels(stimleg{stimg})
    ylabel('ECoG LFB power (x-fold decrease)')
    xlabel(stimgrnames{stimg})

    %exportgraphics(gca,[temp_dir sprintf('task%d.eps',tsk)], 'BackgroundColor','none','ContentType','vector', 'Resolution',300);

end



clear temp*

