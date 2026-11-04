%% Behavior Analysis (Fig. S1) 
% 
% Compute reaction time and accuracy for button presses 
%
% ZQ 2023 Mayo Clinic
%
%%
%% CATEGORIZATION TASK BEHAVIOR ANALYSIS %%
%%
% This is an independent script so setup again
clc; clear; close all

dataPath = setLocalDataPath(1);
subjects = {'sub-01' 'sub-02'};
%
%% Reaction Time Analysis

% total stimuli in the dataset
numstimuli = 24;

% number of trials for each
numreps = 6;

% init varables
reactiontime  = nan(length(subjects),numreps,numstimuli);
resp          = nan(length(subjects),numreps,numstimuli);
buttonpress = zeros(length(subjects),numreps,numstimuli);
stimc       = zeros(length(subjects),numstimuli); 

for zzz = 1%:length(subjects)

    % find behavioral (Log) files
    files = matchfiles( fullfile( dataPath.input, 'derivatives/LogFiles', subjects{zzz}, '*.mat'));
    assert(length(files)==4);

    % process each run (FCFC)
    for ii = 1 : 4
      
      if mod2(ii,2) == 2  % only categoriaztion
      
      % load behavioral file
      a1 = load(files{ii});
      
      % extract the button press time and which button was pressed
      % atlaeast 0.25 s between presses 
      [keytimes,~,keybuttons] = ptviewmoviecheck(a1.timeframes,a1.timekeys,0.25,'t',0.25,1);
      
      for jj = 1 : size(a1.trialpattern, 1)  % 68 trials
            
          ix = find(a1.trialpattern(jj,:));
            
          if ~isempty(ix)  % the first and last are blank
                
              stimix = a1.classorder(ix);  % 1-24 (which stimulus are we on)
              
              fpt = length(a1.timeframes)/size(a1.trialpattern,1);  % number of frames per trial
              assert(fpt==40);
              
              starttime = a1.timeframes((jj-1)*fpt + 1);
              endtime   = a1.timeframes((jj-1)*fpt + fpt + 1);
              
              % keypresses within the stimulus presentation
              okix = find(keytimes > starttime & keytimes < endtime);
              
              stimc(zzz,stimix) = stimc(zzz,stimix) + 1;  % which trial number are we on now
              
              if ~isempty(okix)

                    % mark valid response
                    buttonpress( zzz, stimc(zzz,stimix), stimix) = 1;

                    % reaction time
                    reactiontime(zzz, stimc(zzz,stimix), stimix) = 1000*(keytimes(okix(1)) - starttime);  % just take the first one and covert into ms

                    % which button
                    resp(        zzz, stimc(zzz,stimix), stimix) = keybuttons{okix(1)};
              end
          end
      end
      end
    end
end

clear zzz ii jj a1 ix stimix stimc fpt endtime starttime okix


%% Accuracy

% 1=WORD, 2=FACE, 3=NEITHER

% response key corresponding to the 24 stimuli
resp_key = [0     3     0     1     2     3     3     1     1     3     3     2     2     1     1     1     2     2     2    3     3     3     3     3];

% sub1 uses 55 56 57
% sub2 uses 49 50 51

s1vals = [56 55 57];
s2vals = [49 50 51];

resp_acc = [];

for ii = 1 : length(resp_key)

    % ignoring first and third stim (blank)
    if ii==1 || ii==3
        continue;
    end

    resp_acc(1,ii) = sum( resp(1,:,ii) == s1vals(resp_key(ii)) ) / sum(buttonpress(1,:,ii)) * 100;
    resp_acc(2,ii) = sum( resp(2,:,ii) == s2vals(resp_key(ii)) ) / sum(buttonpress(2,:,ii)) * 100;

end

clear ii
%% Plot reaction time & accuracy

stimgroups  = {[ 6      7      8      9       4  ] [10     11     12     13       5  ] [14    15     16      4   ] [17    18     19      5   ]};
stimleg     = {["0%", "25%", "50%", "75%", "100%"] ["0%", "25%", "50%", "75%", "100%"] ["3%", "5%",  "8%", "100%"] ["4%", "6%", "10%", "100%"]};
stimgrnames = {'WordPhaseCoherence'               'FacePhaseCoherence'               'WordContrast'              'FaceContrast'};

% organize reaction time wrt to stimulus groups
for stimg = 1 : length(stimgroups)

    % temp_dir = fullfile( dataPath.output, 'eps/behavior', stimgrnames{stimg});
    % if not(isfolder(temp_dir))
    %     mkdir(temp_dir)
    % end

    figure

    for zzz = 1:2  % 2 subjects

        temp_yMean =[];
        temp_ySEM  =[];
        temp_acc   =[];
        temp_posx  =[];
        posx = 0;
    
        % for stimuli within the selected stimulus group
        for stimc = 1:length(stimgroups{stimg})

            tempData = squeeze( reactiontime( zzz, :, stimgroups{stimg}(stimc)));
        
            % Mean, SEM and 68% CI
            nTrials = length( tempData) - sum( isnan( tempData));
            yMean = median( tempData, 'omitnan') ;
            ySEM  =  std( tempData, 'omitnan') / sqrt(nTrials);
            
            ts = tinv([0.16, 0.84], nTrials - 1);

            temp_yMean = [temp_yMean yMean];
            temp_ySEM  = [temp_ySEM  (ts(1)*ySEM)'];

            % x-axis locations
            posx = posx + 1;
            temp_posx = [temp_posx posx];

            temp_acc = [temp_acc resp_acc( zzz, stimgroups{stimg}(stimc))];

        end

        % plot reaction time
        subplot (2, 2, 2*zzz-1)
        errorbar(temp_posx, temp_yMean, temp_ySEM, Color = [0 1 0], LineStyle = '-', CapSize=0);

        xlim([0.5 stimc+0.5])
        set(gca,'XTick',1:stimc);
        xticklabels(stimleg{stimg})
        ylabel('Reaction Time (ms)')
        xlabel(stimgrnames{stimg})

        % y-axis range for visualization purposes
        if zzz == 1

            ylim([600 1400])
            title('Sub1')

        elseif zzz == 2

            ylim([1400 2200])
            title('Sub2')

        else

            assert( 1==0, 'Error: only 2 subjects')

        end

        %exportgraphics(gca,[temp_dir sprintf('/RT_sub%d.eps',zzz)], 'BackgroundColor','none','ContentType','vector', 'Resolution',300);

        % plot accuracy
        subplot (2, 2, 2*zzz)
        plot(temp_posx, temp_acc, Color = [0.5 0 0], LineStyle = '-');
        ylabel('Accuracy (%)')
        xlabel(stimgrnames{stimg})

        xlim([0.5 stimc+0.5])
        ylim([50 100])
        set(gca,'XTick',1:stimc);
        xticklabels(stimleg{stimg})
        % show a line for 83.33% accuracy which corresponds to 5 of 6 correct trials
        yline(83.3333, 'k--')

        %exportgraphics(gca,[temp_dir sprintf('/Acc_sub%d.eps',zzz)], 'BackgroundColor','none','ContentType','vector', 'Resolution',300);
    
    end

end

clear temp* zzz posx yMean ySEM ts nTrials

           
%%
%% FIXATION TASK BEHAVIOR ANALYSIS %%
%%
%% Reaction Time Analysis
% similar as in for categorizaiotn

reactiontime  = [];
%resp          = [];
buttonpress = nan(2,154);
stimc       = []; 

for zzz = 1:length(subjects)

    files = matchfiles( fullfile( dataPath.input, 'derivatives/LogFiles', subjects{zzz}, '*.mat'));
    assert(length(files)==4);

    counter = 1;

    % process each run (FCFC)
    for ii = 1 : 4
      
      if mod2(ii,2) == 1  % only fixation
      
      % load behavioral file
      a1 = load(files{ii});
      
      [keytimes,~,keybuttons] = ptviewmoviecheck(a1.timeframes,a1.timekeys,0.25,'t',0.25,1);

      % fixcolor = -a1.fixationorder(find(diff(a1.fixationorder) ~= 0));     % takes the last value before flip

      % a1.fixationorder == -1 ==> red fixation is ON.
      % Since the duration is 600 ms, it's sampled 6 times
      % redonsets are samples when fixation dot is turned red
      redonsets = a1.timeframes( 6 * (find(a1.fixationorder(1:6:end) == -1) - 1) + 1);

      if redonsets(1) == 0
          redonsets = redonsets(2:end);     % ignore the starting 
      end


      % First 8 seconds is rest. Sample rate 2000
      stimonsets = [0 a1.timeframes((4:135)*20 + 1)];
      % Which stimulus (image)?
      stimclass  = [0 flatten([a1.stimclassrec; zeros(size(a1.stimclassrec))])];
 
      for jj = 1 : length(redonsets)  


          % min difference between 2 red onsets is 1.1999

          okix = find(keytimes > redonsets(jj) & keytimes < redonsets(jj)+1.198); 
          if ~isempty(okix)
              % times of more than 1 presses
              % if(length(okix)>1)
              %     ii
              %     length(okix)
              % end
              buttonpress( zzz, counter) = 1;
              reactiontime(zzz, counter) = 1000*(keytimes(okix(1)) - redonsets(jj));  % just take the first one
              %resp(        zzz, counter) = keybuttons{okix(1)};
          else
              buttonpress( zzz, counter) = 0;
              reactiontime(zzz, counter) = nan;
              %resp(        zzz, counter) = 'none';
          end

          stix = find(stimonsets > redonsets(jj));
          if ~isempty(stix)
              stimrec(zzz, counter) = stimclass(stix(1)-1);
          else
              stimrec(zzz, counter) = nan;
          end

          counter = counter+1;
      end
      end
    end
end

clear zzz ii jj a1 okix stix key*


%% Plot  reaction time

stimgroups  = {[ 6      7      8      9       4  ] [10     11     12     13       5  ] [0  14    15     16      4   ] [0  17    18     19      5   ]};
stimleg     = {["0%", "25%", "50%", "75%", "100%"] ["0%", "25%", "50%", "75%", "100%"] ["0%" "3%", "5%",  "8%", "100%"] ["0%", "4%", "6%", "10%", "100%"]};
stimgrnames = {'WordPhaseCoherence'               'FacePhaseCoherence'               'WordContrast'              'FaceContrast'};

% organize reaction time wrt to stimulus groups
for stimg = 1 : length(stimgroups)

    % temp_dir = fullfile( dataPath.output, 'eps/behavior', stimgrnames{stimg});
    % if not(isfolder(temp_dir))
    %     mkdir(temp_dir)
    % end
    
    for zzz = 1:length(subjects)  % 2 subjects

        figure

        temp_yMean =[];
        temp_ySEM  =[];
        temp_posx  =[];
        posx = 0;
    
        % for stimuli within the selected stimulus group
        for stimc = 1:length(stimgroups{stimg})

            ix = find( stimrec(zzz,:) == stimgroups{stimg}(stimc));

            tempData = squeeze( reactiontime( zzz, ix));

            % ignore missed cases
            tempData = tempData(tempData>0);
        
            % Mean, SEM and 68% CI
            nTrials = length( tempData) - sum( isnan( tempData));
            yMean = median( tempData, 'omitnan') ;
            ySEM  =  std( tempData, 'omitnan') / sqrt(nTrials);
            
            ts = tinv([0.16, 0.84], nTrials - 1);

            temp_yMean = [temp_yMean yMean];
            temp_ySEM  = [temp_ySEM  (ts(1)*ySEM)'];

            % x-axis locations
            posx = posx + 1;
            temp_posx = [temp_posx posx];


        end

        % plot reaction time
        errorbar(temp_posx, temp_yMean, temp_ySEM, Color = [0 1 0], LineStyle = '-', CapSize=0);

        xlim([0.5 stimc+0.5])
        set(gca,'XTick',1:stimc);

        % y-axis range for visualization purposes
        if zzz == 1

            ylim([150 500])
            title('Sub1')

        elseif zzz == 2

            ylim([750 1100])
            title('Sub2')

        else

            assert( 1==0, 'Error: only 2 subjects')

        end

        xlabel(stimgrnames{stimg})
        xticklabels(stimleg{stimg})
        ylabel('Reaction Time (ms)')

        %exportgraphics(gca,[temp_dir sprintf('/RTfix_sub%d.eps',zzz)], 'BackgroundColor','none','ContentType','vector', 'Resolution',300);

        
    end

end

clear temp* zzz posx yMean ySEM ts nTrials ix