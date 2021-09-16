%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETUP

%path Documents\git_repo\readingECoG\code

% define
subjects = {'YBA' 'YBD'};
subjectfilenums = {[12 13 14 15] [30 31 32 33]};  % this is for the file names (one file per run per [subject])
channels = 1:128;
numchannels = 128;
photodiode = 129;   % photodiode signal is in #129
epochrng = [-.5 3.5];  % pull out this time range (in seconds)
onsetix = {{1:66 1:66 1:66 1:66} {1:66 1:66 1:66 1:66}};  % indices to pull the black circle (based on the no of stimulus trials per run; photodiode detects this circle)
fsorig = 2000;    % sampling rate (original)
fsjump = 10;      % moving average length
fs = fsorig/fsjump;         		% sampling rate (down-sampled to)
fupper = 200;     % upperbound for frequency analysis
numtasks = 2;     % alternate the tasks
numreps = 6;      % 6 total trials for each numtasks (FC#1: 3 trials, FC#2: 3 trials)
numstimuli = 24;  % remember that we don't use #1 and #3, so entries for this will be blank 
numtrials = 66;   % in one run, there are these many stimulus trials (numstimuli (=22) X 3 trials per numstimuli)
numruns = 4;      % total number of runs per subject (FCFC)
tasklabels = {'Fixation' 'Categorization'};
 
% calc
epochtime_bb = (epochrng(1)*fs : epochrng(2)*fs)/fs;
epochtime = (epochrng(1)*fsorig : epochrng(2)*fsorig)/fsorig;  % time in seconds for each data point of an epoch
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETUP CONTINUED (ASSIGN LABELS TO CHANNELS)

% define
cfun = @(x,y) cellfun(@(z) [x z],mat2cellstr(y),'UniformOutput',0); % mat2cellstr: return a cell matrix of strings (knkutils)

% initialize
channellabels = {};  % 128 channels for each subject!

% YBA:
channellabels{2} = [ ...
cfun('LAT',1:4) ...
cfun('LAIT',1:4) ...
cfun('LPIT',1:4) ...
    cfun('LPT',1:6) ...
cfun('LD',1:10) ...
cfun('RF',1:6) ...
cfun('RAIT',1:4) ...
cfun('RMIT',1:4) ...
cfun('RPIT',1:4) ...
cfun('RATD',1:10) ...
cfun('RPTD',1:8) ...
cfun('RP',1:6) ...
cfun('RPTO',1:8) ...
cfun('RMini',1:16) ...
cfun('EKG',1:2) ...
repmat({'nolabel'},[1 32]) ...
];

% initialize
anatlabels = {};  % 128 channels for each subject!

% YBA:
anatlabels{2} = [ ...
repmat({'AT'},[1 4]) ...   % AT means anterior temporal (on the ventral aspect)
repmat({'AT'},[1 4]) ...
repmat({'IT'},[1 4]) ...   % IT means inferior temporal (on the ventral aspect)
repmat({'LL'},[1 6]) ...   % LL means lateral lateral stuff. like STS stuff
repmat({'DD'},[1 10]) ...  % DD means weird depth electrodes (hippocampal?)
repmat({'F'}, [1 6]) ...   % F means frontal
repmat({'AT'},[1 4]) ...
repmat({'IT'},[1 4]) ...
repmat({'IT'},[1 4]) ...
repmat({'DD'},[1 10]) ...
repmat({'DD'},[1 8]) ...
repmat({'P'},[1 6]) ...    % P means parietal
repmat({'O'},[1 8]) ...    % O means occipital (could be early visual or lateral visual)
repmat({'M'},[1 16]) ...   % M means mini (i think this also means they are "O", occipital)
repmat({'X'},[1 2]) ...    % X means junk (EKG or nolabel)
repmat({'X'},[1 32]) ...
];

%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREP DATA

% do it

data = zeros(length(epochtime),length(subjects),numchannels,numtasks,numreps,numstimuli,'single');
bb_data = zeros(length(epochtime_bb),length(subjects),numchannels,numtasks,numreps,numstimuli,'single');
nb_data = zeros(length(epochtime_bb),length(subjects),numchannels,numtasks,numreps,numstimuli,'single');
alpha_data = zeros(length(epochtime_bb),length(subjects),numchannels,numtasks,numreps,numstimuli,'single');
beta_data = zeros(length(epochtime_bb),length(subjects),numchannels,numtasks,numreps,numstimuli,'single');
theta_data = zeros(length(epochtime_bb),length(subjects),numchannels,numtasks,numreps,numstimuli,'single');
delta_data = zeros(length(epochtime_bb),length(subjects),numchannels,numtasks,numreps,numstimuli,'single');
psd_on = zeros(length(subjects),numchannels,numtasks,numreps,numstimuli,fupper,'single');
psd_off = zeros(length(subjects),numchannels,numtasks,numreps,numstimuli,fupper,'single');
stimcounter = zeros(length(subjects),numchannels,numtasks,numstimuli);  % number of trials encountered so far
recdata = {};
bbtemp = {};


%% Finding epochs

onsets = zeros(numtrials,length(subjects),numruns);  % indices of data points at which trials occur
onsetpoints = {};


for zzz=1:length(subjects)

  % get behavioral files
  files = matchfiles(sprintf('../data/%s/LogFiles/*.mat',subjects{zzz}))   %matchfiles: to match the filename (knkutils)
  assert(length(files)==length(subjectfilenums{zzz}));
  
  temponsetpoints = []; 
  optemp = 0;
  % process each run (photodiode)
  for p=1:length(subjectfilenums{zzz})

    % load photodiode information
    chantoload = photodiode;     % the special 129
    file0 = sprintf('../data/%s/%sDatafile%03d_ch%d.mat', ...
                    subjects{zzz},subjects{zzz},subjectfilenums{zzz}(p),chantoload);
    pd = load(file0);
    fs0 = pd.analogInfos.SampleRate; assert(fs0==fsorig);
    
    % get good time region (rough pass)
    numsam = length(pd.analogTraces);
    [~,theseOnsets] = findpeaks(-pd.analogTraces,'MinPeakDistance',3*fsorig,'MinPeakHeight',150);  % if black is low, then - means find pos deflections
    theseOnsets = theseOnsets(onsetix{zzz}{p});  % this is not quite right (finding peaks), but roughly right
    
    % extract good region
    okixs = theseOnsets(1)-fsorig*1 : theseOnsets(end)+fsorig*1;  % 1 s padding on both sides
    tempdata = pd.analogTraces(okixs);
    tempdata = tsfilter(tempdata,constructbutterfilter1D(length(tempdata),-round(numtrials/2)));  % freq higher than N cycles per run are passed  %%constructbutterfilter1D, tsfilter (knkutils)
    locoftransition = find(diff(sign(tempdata))==-2);  % this assigns to the leftward side of the transition... oh well. small (0.5 ms).
    theseOnsets = (okixs(1) - 1) + locoftransition;
    assert(length(theseOnsets)==numtrials);  % sanity check that we got the right number of stimulus trials

    % visualize for sanity
    %figureprep([100 100 1000 300]); hold on;		% figureprep & figurewrite(knkutils)
    figureprep([100 100 1000 300]); hold on;
    plot(pd.analogTraces);
    straightline(theseOnsets,'v','m-');
    figurewrite(sprintf('photodiode_subj%d_file%d',zzz,p),[],[],'inout/photodiode');

    % record
    onsets(:,zzz,p) = theseOnsets;
    temponsetpoints = [temponsetpoints optemp+theseOnsets]; % cumulative epoch points
    optemp = optemp + length(pd.analogTraces);	% cumulative length of runs 
  end
  onsetpoints{zzz} = single(temponsetpoints);
  

end

%%  Channel Analysis


for zzz=1:length(subjects)
    
  % get behavioral files
  files = matchfiles(sprintf('../data/%s/LogFiles/*.mat',subjects{zzz}))   %matchfiles: to match the filename (knkutils)
  assert(length(files)==length(subjectfilenums{zzz}));
  
  % process each channel
  for ccc=1:numchannels
    % init
    
    sprintf('Subject%d_%03d',zzz,ccc)
    collectdata = [];
    bblengths = [];
    
    % process each run (actual data)
    for p=1:length(subjectfilenums{zzz})

      % load data
      chantoload = channels(ccc);  % the usual 1-128
      file0 = sprintf('../data/%s/%sDatafile%03d_ch%d.mat', ...
                      subjects{zzz},subjects{zzz},subjectfilenums{zzz}(p),chantoload);
      if ~exist(file0,'file')
        continue;
      end
      pd = load(file0);
      fs0 = pd.analogInfos.SampleRate; assert(fs0==fsorig);
      
      % record data
      collectdata = [collectdata pd.analogTraces];
      bblengths(p) = length(pd.analogTraces);

    end
    
    % if no data exists, get out early
    if isempty(collectdata)
      continue;
    end
    
    % now do broadband:
    
	% OLD way (version 1): bb = g(pd.analogTraces);   % <===== NOTE THIS!!! 

    collectdata = collectdata';  %convert to time X channels
   
    % TOTAL HACK SPIKE (detect later!)
    if zzz==1   %subject YBA
      badrng = 1.072 * 10^6 : 1.096 * 10^6; % data seems to be corrupted in this time range for all channels
      %bb(badrng) = NaN;  %%median(bb(setdiff(1:length(bb),badrng)));
      collectdata(badrng) = NaN;
      %collectdata(badrng) = nanmedian(collectdata);
    end
    
	% NEW way:
    bands = [70  90;  % 20 hz bins, avoiding 60 and 180 (max 210)
             90  110;
             110 130;
             130 150;
             150 170];   %% HACK OUT
%             190 210];
    narrowband = [30 50];  
    beta  = [12 30];
    alpha = [8 12];
    theta = [4 7];
    delta = [1 4];
    
    %check the bb influence on nb
    
    bb = ecog_extractBroadband(collectdata,fsorig,[],bands);  % NOTE: WE DO WHOLE EXPERIMENT AT ONCE TO ENSURE EQUAL SCALING ISSUE...
                                                              % ecog_extractBroadband: mean after hilbert; power; geomean (ECoG utilities by JW_NYU)
    
    nb = ecog_extractBroadband(collectdata,fsorig,[],narrowband); 
    betab  = ecog_extractBroadband(collectdata,fsorig,[],beta);
    alphab = ecog_extractBroadband(collectdata,fsorig,[],alpha);
    thetab  = ecog_extractBroadband(collectdata,fsorig,[],theta);
    deltab = ecog_extractBroadband(collectdata,fsorig,[],delta);
    
    %bb_bp = ieeg_butterpass(collectdata, [70 170], fsorig);   % bandpass; amplitude (mnl_ieegBasics by DH)

	% 	RAW CASE: bb = pd.analogTraces;

    
    bbtemp{zzz,ccc} = single(bb);   % save broadband so we can inspect it!
    recdata{zzz,ccc} = single(collectdata);    % raw data amplitude

    
    %assert(all(isfinite(bb)));
    
    % compute moving average to reduce data rate to 200Hz (2000/10)
    bb = movmean(bb,  fsjump,1);  % slight slop due to evenness of the number.  we can plop this in the photodiode transition issue.
    nb = movmean(nb,  fsjump,1);
    betab = movmean(betab,  fsjump,1);
    alphab = movmean(alphab,  fsjump,1);
    thetab = movmean(thetab,  fsjump,1);
    deltab = movmean(deltab,  fsjump,1);    
    
    % process each run
    for p=1:length(subjectfilenums{zzz})

      % load behavioral file
      a1 = load(files{p});

      % extract epochs (trials)
      for ttt=1:size(onsets,1)
          
        % indices 
        ix_bb = onsets(ttt,zzz,p) + (epochrng(1)*fsorig : fsjump : epochrng(2)*fsorig);  % NOTE: the 10 makes things nice and even, and hits 0
        ix_raw = onsets(ttt,zzz,p) + (epochrng(1)*fsorig : epochrng(2)*fsorig);	% channel analysis	 (t = -0.5:3.5 s; f = 2 kHz)
        on_samples  = length(epochrng(1)*fsorig:0*fsorig) : length(epochrng(1)*fsorig:0*fsorig) + 2 * fsorig - 1;  % samples corresponding to on stimuli (0 - 2 s)
        %off_samples = cat (2, 1 : length(epochrng(1)*fsorig:0*fsorig)-1, length(epochrng(1)*fsorig:0*fsorig) + 2 * fsorig : 4*fsorig);  
        off_samples = 1 : length(epochrng(1)*fsorig:0*fsorig)-1;  % samples corresponding to off-stimuli-period (-0.5 - 0 s)
        
        temp_bb = bb(sum(bblengths(1:p-1)) + ix_bb);			% bb analysis
        temp_nb = nb(sum(bblengths(1:p-1)) + ix_bb);            % nb analysis
        temp_beta = betab(sum(bblengths(1:p-1)) + ix_bb);       % beta band analysis
        temp_alpha = alphab(sum(bblengths(1:p-1)) + ix_bb);     % alpha band analysis
        temp_theta = thetab(sum(bblengths(1:p-1)) + ix_bb);     % alpha band analysis
        temp_delta = deltab(sum(bblengths(1:p-1)) + ix_bb);     % alpha band analysis
        
		temp_raw = collectdata(sum(bblengths(1:p-1)) + ix_raw);	% raw analysis
%         if any(isnan(temp_bb))
%           fprintf('BADDATA: ttt=%d, p=%d, ccc=%d, zzz=%d\n',ttt,p,ccc,zzz);
%           temp_bb(:) = badval;
%         end
%         if any(isnan(temp_raw))
% 		  temp_raw(:) = badval;
%         end
        
		stimco = stimcounter(zzz,ccc,mod2(p,2),a1.stimclassrec(ttt))+1;   % a1.stimclassrec tells us the stim number (1-24)
        
        bb_data(:,zzz,ccc,mod2(p,2), stimco, ...
             a1.stimclassrec(ttt)) = temp_bb;
        nb_data(:,zzz,ccc,mod2(p,2), stimco, ...
             a1.stimclassrec(ttt)) = temp_nb;
        beta_data(:,zzz,ccc,mod2(p,2), stimco, ...
             a1.stimclassrec(ttt)) = temp_beta;
        alpha_data(:,zzz,ccc,mod2(p,2), stimco, ...
             a1.stimclassrec(ttt)) = temp_alpha;
        theta_data(:,zzz,ccc,mod2(p,2), stimco, ...
             a1.stimclassrec(ttt)) = temp_theta;
        delta_data(:,zzz,ccc,mod2(p,2), stimco, ...
             a1.stimclassrec(ttt)) = temp_delta;
         
         
        data(:,zzz,ccc,mod2(p,2), stimco, ...
             a1.stimclassrec(ttt)) = temp_raw;
        
%         [psdvar,f] = pwelch(temp_raw(on_samples)',hamming(1000),0,2 ^ nextpow2(fsorig),fsorig);
%          %figureprep([100 100 900 300],1);plot(f,10*log10(psdvar));
%         psd_on(zzz,ccc,mod2(p,2), stimco, ...
%              a1.stimclassrec(ttt),:) = psdvar(1:fupper);
%         [psdvar,f] = pwelch(temp_raw(off_samples)',hamming(1000),0,2 ^ nextpow2(fsorig),fsorig);
%         psd_off(zzz,ccc,mod2(p,2), stimco, ...
%              a1.stimclassrec(ttt),:) = psdvar(1:fupper);
		
        stimcounter(zzz,ccc,mod2(p,2),a1.stimclassrec(ttt)) = ...
           stimcounter(zzz,ccc,mod2(p,2),a1.stimclassrec(ttt)) + 1;
      end
    
    end
  end  
end

%% CHANNEL AVG TIME AND FREQUENCY RESPONSE

psd_avg = squeeze(mean(reshape(psd_on,numchannels*length(subjects),numtasks*numreps*numstimuli,fupper),2));
psd_baseline = squeeze(mean(reshape(psd_off,numchannels*length(subjects),numtasks*numreps*numstimuli,fupper),2));

bb_data_avg = mean(reshape(bb_data,length(epochtime_bb),length(subjects)*numchannels,numtasks*numreps*numstimuli),3)';
nb_data_avg = mean(reshape(nb_data,length(epochtime_bb),length(subjects)*numchannels,numtasks*numreps*numstimuli),3)';
alpha_data_avg = mean(reshape(alpha_data,length(epochtime_bb),length(subjects)*numchannels,numtasks*numreps*numstimuli),3)';
beta_data_avg = mean(reshape(beta_data,length(epochtime_bb),length(subjects)*numchannels,numtasks*numreps*numstimuli),3)';
data_avg = mean(reshape(data,length(epochtime),length(subjects)*numchannels,numtasks*numreps*numstimuli),3)';

%% PLOT CHANNEL AVG TIME AND FREQUENCY RESPONSE

clims = [0 50];
figure,imagesc(bb_data_avg);
colorbar;
set(gca,'XTick',0:100:800)
set(gca,'XTickLabel',-0.5:0.5:3.5)

% EP odd-sub1/even-sub2
figure,imagesc(data_avg,[-50,50]);
colorbar;
set(gca,'XTick',0:1000:8000)
set(gca,'XTickLabel',-0.5:0.5:3.5)

%%
for ccc = 41
    figure,plot(data_avg(ccc*2-1,:))
end
hold on;
plot(data_avg(76*2,:))
plot(data_avg(77*2,:))
plot(data_avg(78*2,:))
hold off;
set(gca,'XTick',0:1000:8000)
set(gca,'XTickLabel',-0.5:0.5:3.5)

%% Plot normalized PSD

for ccc = 41
    figure,semilogy(psd_avg(ccc*2-1,:),'r')
    hold on;
    semilogy(psd_baseline(ccc*2-1,:),'k')
    hold off;
end
plot((normalize(10*log10(psd_avg(ccc,:)) - 10*log10(psd_baseline(ccc,:)))))

%% Good Channels(Visually Identified)

%EVCgcc = {horzcat(74:75),horzcat(43:44, 108)};

%LVCgcc = {horzcat(76:78),horzcat(45:46, 109:110)};

%gcc = sort(horzcat((9:11)*2-1, (44:46)*2-1, (74:78)*2-1, (3)*2, (20:22)*2, (86:87)*2, (108:110)*2));

gcc = {horzcat((9:11), (44:46), (74:78)), horzcat((3), (20:22), (86:87), (108:110)) };
%% Stimuli Timecourse

stimgroups  = {[6 7 8 9 4]   [10 11 12 13 5] [14 15 16 4] [17 18 19 5]};% [20 21 22 10] [2 23 24]};
stimleg  = {["0", "25", "50", "75", "100"]  ["0", "25", "50", "75", "100"] ["3" "5" "8" "100"] ["4" "6" "10" "100"]};
stimgrnames = {'Word Phase' 'Face Phase'    'Word Con'   'Face Con'};%   'Noise Con'   'Other'};
stimresp_f = zeros(length(epochtime_bb),length(subjects),20,length(stimgroups),5);
stimresp_c = zeros(length(epochtime_bb),length(subjects),20,length(stimgroups),5);
stimresp_g = zeros(length(subjects),20,2,length(stimgroups),5);
counter = 1;

for zzz = 1:length(subjects)
   
    for ccc = gcc{zzz}
        
        figureprep([100 100 1700 1100]);
    
        for stimg = 1:length(stimgroups)
            
            temp = [];
            subplot(4,length(stimgroups),stimg)
            for stimc = 1:length(stimgroups{stimg})
                
%                 temp = [temp smooth(squeeze(mean(bb_data(:,zzz,ccc,1,:,stimgroups{stimg}(stimc)),5)),40)];

                stimresp_f(:,zzz,counter,stimg,stimc) = smooth(squeeze(mean(bbdata_br(:,zzz,ccc,1,:,stimgroups{stimg}(stimc)),5)),40);
                
                plot(smooth(squeeze(mean(bbdata_br(:,zzz,ccc,1,:,stimgroups{stimg}(stimc)),5)),40));
                ylim([-50 200]); 
                xlim([0 400]);
                set(gca,'XTick',0:100:400); set(gca,'XTickLabel',-0.5:0.5:1.5);
                hold on;
                
            end
            %stimresp_f{stimg,stimc} = temp;    
            hold off;
            title(sprintf('%s_F',stimgrnames{stimg}));
            legend(stimleg{stimg});
            
            temp = [];
            subplot(4,length(stimgroups),length(stimgroups)+stimg)
            for stimc = 1:length(stimgroups{stimg})
                
                %temp = [temp smooth(squeeze(mean(bb_data(:,zzz,ccc,2,:,stimgroups{stimg}(stimc)),5)),40)];
                stimresp_c(:,zzz,counter,stimg,stimc) = smooth(squeeze(mean(bbdata_br(:,zzz,ccc,2,:,stimgroups{stimg}(stimc)),5)),40);
                plot(smooth(squeeze(mean(bbdata_br(:,zzz,ccc,2,:,stimgroups{stimg}(stimc)),5)),40));
                ylim([-50 200]); 
                xlim([0 400]);
                set(gca,'XTick',0:100:400); set(gca,'XTickLabel',-0.5:0.5:1.5);
                hold on;
                
            end
            %stimresp_c{stimg} = temp;
            hold off;
            title(sprintf('%s_C',stimgrnames{stimg}));
            legend(stimleg{stimg});
            
            temp = [];
            subplot(4,length(stimgroups),2*length(stimgroups)+stimg)
            for stimc = 1:length(stimgroups{stimg})-1
                stimresp_g(zzz,ccc,1,stimg,stimc) = mean(mean(bbdata_br(101:300,zzz,ccc,1,:,stimgroups{stimg}(stimc)),5),1);
                plot([stimc,stimc+1],[mean(mean(bbdata_br(101:300,zzz,ccc,1,:,stimgroups{stimg}(stimc)),5),1),mean(mean(bbdata_br(101:300,zzz,ccc,1,:,stimgroups{stimg}(stimc+1)),5),1)],'r-', 'LineWidth', 2);
                ylim([-50 100]);
                hold on;
                stimresp_g(zzz,ccc,2,stimg,stimc) = mean(mean(bbdata_br(101:300,zzz,ccc,2,:,stimgroups{stimg}(stimc)),5),1);
                plot([stimc,stimc+1],[mean(mean(bbdata_br(101:300,zzz,ccc,2,:,stimgroups{stimg}(stimc)),5),1),mean(mean(bbdata_br(101:300,zzz,ccc,2,:,stimgroups{stimg}(stimc+1)),5),1)],'b-', 'LineWidth', 2);
                hold on;
            end
            stimc = stimc + 1;
            stimresp_g(zzz,ccc,1,stimg,stimc) = mean(mean(bbdata_br(101:300,zzz,ccc,1,:,stimgroups{stimg}(stimc)),5),1);
            stimresp_g(zzz,ccc,2,stimg,stimc) = mean(mean(bbdata_br(101:300,zzz,ccc,2,:,stimgroups{stimg}(stimc)),5),1);
            %stimresp_c{stimg} = temp;
            hold off;
            title(sprintf('%s_G',stimgrnames{stimg}));
            legend('F','C')
            
            temp = [];
            subplot(4,length(stimgroups),3*length(stimgroups)+stimg)
            for stimc = 1:length(stimgroups{stimg})
                bar(stimc,nanmean(reactiontime(zzz,:,stimgroups{stimg}(stimc)),2));
                hold on;
            end
            hold off;
            title(sprintf('%s_t_reaction',stimgrnames{stimg}));
            legend(stimleg{stimg});
            
        end
        figurewrite(sprintf('Subj%d_ch%03d',zzz,ccc),[],[],'stimtimecourse');
        
        counter = counter+1;
        
    end
    
end

%%


for stimg = 1:length(stimgroups)
    
    figureprep([100 100 1700 1100]);
    
    for stimc = 1:length(stimgroups{stimg})
        
        subplot(3,length(stimgroups{stimg}),stimc);
        plot(mean(mean(stimresp_f(:,1,7:11,stimg,stimc),2),3));
        ylim([0 50]); 
        xlim([0 400]);
        set(gca,'XTick',0:100:400); set(gca,'XTickLabel',-0.5:0.5:1.5);
        hold on;
        plot(mean(mean(stimresp_c(:,1,7:11,stimg,stimc),2),3));
        hold off;
        legend('F','C');
        
        subplot(4,length(stimgroups{stimg}),length(stimgroups{stimg})+stimc);
        plot(mean(mean(stimresp_c(:,1,7:11,stimg,stimc),2),3)-mean(mean(stimresp_f(:,1,7:11,stimg,stimc),2),3));
        ylim([0 50]); 
        xlim([0 400]);
        [h,p,ci,stats] = ttest2(mean(mean(stimresp_c(:,1,7:11,stimg,stimc),2),3),mean(mean(stimresp_f(:,1,7:11,stimg,stimc),2),3));
        title(num2str([h,p]));
        set(gca,'XTick',0:100:400); set(gca,'XTickLabel',-0.5:0.5:1.5);
    end
    
    subplot(3,length(stimgroups{stimg}),2*length(stimgroups{stimg})+1);
    
    for stimc = 1:length(stimgroups{stimg})-1
        
        plot([stimc,stimc+1],[mean(mean(stimresp_g(:,:,1,stimg,stimc),2),1),mean(mean(stimresp_g(:,:,1,stimg,stimc+1),2),1)],'r-', 'LineWidth', 2);
        ylim([-50 100]);
        hold on;
        plot([stimc,stimc+1],[mean(mean(stimresp_g(:,:,2,stimg,stimc),2),1),mean(mean(stimresp_g(:,:,2,stimg,stimc+1),2),1)],'b-', 'LineWidth', 2);
        hold on;
    end
    hold off;
    legend('F','C')
    
    subplot(3,length(stimgroups),2*length(stimgroups{stimg})+2)
    for stimc = 1:length(stimgroups{stimg})
    bar(stimc,nanmean(nanmean(reactiontime(zzz,:,stimgroups{stimg}(stimc)),2),1));
    hold on;
    end
    hold off;
    title(sprintf('%s',stimgrnames{stimg}));
    figurewrite(sprintf('%s',stimgrnames{stimg}),[],[],'FvsC');
end
%%

    


%%
% what we did:
% - did a moving average of size 10 (to go from 2000 to 200)
% - for each of the 66 trials, define time=0 based on the photodiode zero-crossing (leftward).
%   then extract the epoch range [-.5 3.5]*fs from the broadband time-series and save it.
%   note that these epochs no longer overlap a bit!
%
% outputs:
% - <data> is epochtime x 3 subjects x 128 channels x 2 tasks x 6 trials x 24 stimuli
%   - note that the chronological order of the trials is preserved.
%   - note that stimulus #1 and #3 are not presented, so the data for these will be zeros.
%
% NOTICE HOW WE COMPUTE BB ON THE ENTIRE DATASET AT ONCE!  otherwise, you'll get weird run to run differences.

%%

bbdata_br = bb_data - nanmean(nanmean(bb_data(1:100, :, :, :, :, :),5),1);
bbdata_pc = bsxfun(@rdivide,bb_data,nanmean(nanmean(bbdata_br(:, :, :, :, :, setdiff(1:24,[1 3])),6),4))*100;

% baseline subtraction and normalization (per frequency)

% for zzz = 
%     for ccc = 
%         for ttt =
%             for stim
    [S, f] = getWaveletSpectrogram(squeeze(mean(mean(data(:, 1, 75, 2, :, :),5),6)), fsorig, [1, 200]);     % Returns the Morlet (Gabor) wavelet transform (Spectrogram) for a signal - HH
    %[S2, f] = getWaveletSpectrogram(squeeze(mean(mean(data(off_samples, 1, 75, 2, :, :),5),6)), fsorig, [1, 200]);   
    figure,uimagesc(epochtime,f,S)
    axis xy
    plot(S)

%% fmri type result - beta gain ....  Already included

ccc = 75
for j=1:4
    subplot(1,4,j);
for stimc = 1:length(stimgroups{j})-1
    plot([stimc,stimc+1],[mean(mean(bbdata_br(101:300,1,ccc,1,:,stimgroups{j}(stimc)),5),1),mean(mean(bbdata_br(101:300,1,ccc,1,:,stimgroups{j}(stimc+1)),5),1)],'r-', 'LineWidth', 2);
    ylim([-50 200]);
    hold on;
    plot([stimc,stimc+1],[mean(mean(bbdata_br(101:300,1,ccc,2,:,stimgroups{j}(stimc)),5),1),mean(mean(bbdata_br(101:300,1,ccc,2,:,stimgroups{j}(stimc+1)),5),1)],'b-', 'LineWidth', 2);
    hold on;
end
hold off;
end
%% Reaction Time

reactiontime = NaN * ones(length(subjects),numreps,numstimuli);
stimc = zeros(length(subjects),numstimuli); 
for zzz = 1:length(subjects)
    files = matchfiles(sprintf('../data/%s/LogFiles/*.mat',subjects{zzz}))   %matchfiles: to match the filename (knkutils)
    assert(length(files)==length(subjectfilenums{zzz}));
    % process each run
    for p=1:length(subjectfilenums{zzz})
      if mod2(p,2)==2
      % load behavioral file
      a1 = load(files{p});
      [keytimes,badtimes,keybuttons] = ptviewmoviecheck(a1.timeframes,a1.timekeys,0.25,'t',0.25,1);
          for q=1:size(a1.trialpattern)  % 68 trials
            ix = find(a1.trialpattern(q,:));
            if ~isempty(ix)  % the first and last are blank
                stimix = a1.classorder(ix);  % 1-24 (which stimulus are we on)
                fpt = length(a1.timeframes)/size(a1.trialpattern,1);  % number of frames per trial
                assert(fpt==40);
                starttime = a1.timeframes((q-1)*fpt + 1);
                endtime   = a1.timeframes((q-1)*fpt + fpt + 1);
                okix = find(keytimes > starttime & keytimes < endtime);
                stimc(zzz,stimix) = stimc(zzz,stimix) + 1;  % which trial number are we on now?
                if ~isempty(okix)
                    reactiontime(zzz,stimc(zzz,stimix),stimix) = 1000*(keytimes(okix(1)) - starttime);  % just take the first one
                end
            end
          end  
      end
    end
end




