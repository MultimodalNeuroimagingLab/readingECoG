%% Stimulus and task effects in alpha/beta (Fig.S9) 
% 
% time-average alpha/beta power changes during fixaiton vs categorizaiton
%
% ZQ 2023 Mayo Clinic
%
%%  Alpha - Baseline subtraction and normalization
% similar to HFB (above)

% Removing bad trials (4 in total for subject 1, 3 of 4 during categorization)

% WC 5 - trial 2
alpha_data(:, 1, :, 2, 2, 15) = nan;

% WPC 0 - trial 3
alpha_data(:, 1, :, 2, 3, 6) = nan;

% WPC 50 - trial 3
alpha_data(:, 1, :, 2, 3, 8) = nan;

% WC/WPC 100 - trial 5 - fix
alpha_data(:, 1, :, 1, 5, 4) = nan;

winlen_base = find(epochtime_bb >= -0.5 & epochtime_bb < -0.1);
temp_base = reshape( alpha_data( winlen_base, :, :, :, :, :), length(winlen_base), size(alpha_data,2), size(alpha_data,3), size(alpha_data,4) * size(alpha_data,5) * size(alpha_data,6));
alpha_base = geomean( geomean( temp_base, 4, 'omitnan'), 1, 'omitnan');

% lfb x-fold change
alphadata_pc = bsxfun( @rdivide, alpha_data, alpha_base) - 1;

%%  Beta - Baseline subtraction and normalization
% similar to HFB (above)

% Removing bad trials (4 in total for subject 1, 3 of 4 during categorization)

% WC 5 - trial 2
beta_data(:, 1, :, 2, 2, 15) = nan;

% WPC 0 - trial 3
beta_data(:, 1, :, 2, 3, 6) = nan;

% WPC 50 - trial 3
beta_data(:, 1, :, 2, 3, 8) = nan;

% WC/WPC 100 - trial 5 - fix
beta_data(:, 1, :, 1, 5, 4) = nan;

winlen_base = find(epochtime_bb >= -0.5 & epochtime_bb < -0.1);
temp_base = reshape( beta_data( winlen_base, :, :, :, :, :), length(winlen_base), size(beta_data,2), size(beta_data,3), size(beta_data,4) * size(beta_data,5) * size(beta_data,6));
beta_base = geomean( geomean( temp_base, 4, 'omitnan'), 1, 'omitnan');

% lfb x-fold change
betadata_pc = bsxfun( @rdivide, beta_data, beta_base) - 1;


%%
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


%%
%%  Alpha averaged over time window
% similar to HFB but with slight offset and reverse y-axis
%%

% temp_dir = fullfile( dataPath.output, 'eps/avg/al', roiName);
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
                    squeeze( geomean( alphadata_pc( winlen_resp, zzz, ccc, :, :, stimgroups{stimg}(stimc)) + 1, 1, 'omitnan'))';
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
    ylabel('ECoG Alpha power (x-fold decrease)')
    xlabel(stimgrnames{stimg})

    %exportgraphics(gca,[temp_dir sprintf('/avg%s.eps',stimgrnames{stimg})], 'BackgroundColor','none','ContentType','vector', 'Resolution',300);

end

clear temp*


%%
%%  Beta averaged over time window
% similar to alpha
%%
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
                    squeeze( geomean( betadata_pc( winlen_resp, zzz, ccc, :, :, stimgroups{stimg}(stimc)) + 1, 1, 'omitnan'))';
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
    ylabel('ECoG Beta power (x-fold decrease)')
    xlabel(stimgrnames{stimg})

    %exportgraphics(gca,[temp_dir sprintf('/avg%s.eps',stimgrnames{stimg})], 'BackgroundColor','none','ContentType','vector', 'Resolution',300);

end



clear temp*