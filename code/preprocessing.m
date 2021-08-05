%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETUP


% define
subjects = {'YBA'};
subjectfilenums = {[12 13 14 15]};  % this is for the file names (one file per run per [subject])
channels = 1:128;
numchannels = 96;   % only for subject YBA; when including other subjects make 128;
photodiode = 129;   % photodiode signal is in #129
epochrng = [-.5 3.5];  % pull out this time range (in seconds)
onsetix = {{1:66 1:66 1:66 1:66}};  % indices to pull the black circle (based on the no of stimulus trials per run; photodiode detects this circle)
fsorig = 2000;    % sampling rate (original)
fsjump = 10;      % moving average length
fs = fsorig/fsjump;         		% sampling rate (down-sampled to)
fupper = 200;     % upperbound for frequency analysis
numtasks = 2;     % alternate the tasks
numreps = 6;      % 6 total trials for each numtasks (FC#1: 3 trials, FC#2: 3 trials)
numstimuli = 24;  % remember that we don't use #1 and #3, so entries for this will be blank 
numtrials = 66;   % in one run, there are these many stimulus trials (numstimuli X 3 trials per numstimuli)
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
onsets = zeros(numtrials,length(subjects),numruns);  % indices of data points at which trials occur
data = zeros(length(epochtime),length(subjects),numchannels,numtasks,numreps,numstimuli,'single');
bb_data = zeros(length(epochtime_bb),length(subjects),numchannels,numtasks,numreps,numstimuli,'single');
psd_on = zeros(length(subjects),numchannels,numtasks,numreps,numstimuli,fupper,'single');
psd_off = zeros(length(subjects),numchannels,numtasks,numreps,numstimuli,fupper,'single');
stimcounter = zeros(length(subjects),numchannels,numtasks,numstimuli);  % number of trials encountered so far
recdata = {};
bbtemp = {};
bb_bptemp = {};
for zzz=1:length(subjects)

  % get behavioral files
  files = matchfiles(sprintf('../data/%s/LogFiles/*.mat',subjects{zzz}))   %matchfiles: to match the filename (knkutils)
  assert(length(files)==length(subjectfilenums{zzz}));
  
  onsetpoints = []; % done only for YBA. need to revisit while including other subjects
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
    figure('position',[100 100 1000 300]); hold on;
    plot(pd.analogTraces);
    straightline(theseOnsets,'v','m-');
    figurewrite(sprintf('photodiode_subj%d_file%d',zzz,p),[],[],'~/inout/photodiode');

    % record
    onsets(:,zzz,p) = theseOnsets;
    onsetpoints = [onsetpoints optemp+theseOnsets]; % cumulative epoch points
    optemp = optemp + length(pd.analogTraces);	% cumulative length of runs 
  end
  
  % process each channel
  for ccc=1:numchannels
    % init
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
	% NEW way:
    bands = [70 90;  % 20 hz bins, avoiding 60 and 180 (max 210)
             90 110;
             110 130;
             130 150;
             150 170];   %% HACK OUT
%             190 210];
    bb = ecog_extractBroadband(collectdata,fsorig,[],bands);  % NOTE: WE DO WHOLE EXPERIMENT AT ONCE TO ENSURE EQUAL SCALING ISSUE...
                                                              % ecog_extractBroadband: mean after hilbert; power; geomean (ECoG utilities by JW_NYU)
 
    bb_bp = ieeg_butterpass(collectdata, [70 170], fsorig);   % bandpass; amplitude (mnl_ieegBasics by DH)

	% 	RAW CASE: bb = pd.analogTraces;

    
    bbtemp{zzz,ccc} = single(bb);   % save broadband so we can inspect it!
    bb_bptemp{zzz,ccc} = single(bb_bp);
    

    % TOTAL HACK SPIKE (detect later!)
    assert(all(isfinite(bb)));
    if zzz==1   %subject YBA
      badrng = 1.072 * 10^6 : 1.09 * 10^6; % data seems to be corrupted in this time range for all channels
      bb(badrng) = NaN;  %%median(bb(setdiff(1:length(bb),badrng)));
      collectdata(badrng) = NaN;
      badval = nanmedian(bb);  % **ROBUST
      badvalraw = nanmedian(collectdata);
    end
    
    recdata{zzz,ccc} = single(collectdata);    % raw data amplitude
    % compute moving average to reduce data rate to 200Hz (2000/10)
    bb = movmean(bb,  fsjump,1);  % slight slop due to evenness of the number.  we can plop this in the photodiode transition issue.
    bb_bp = movmean(bb_bp,  fsjump,1);
    
    % process each run
    for p=1:length(subjectfilenums{zzz})

      % load behavioral file
      a1 = load(files{p});

      % extract epochs (trials)
      for ttt=1:size(onsets,1)
        ix_bb = onsets(ttt,zzz,p) + (epochrng(1)*fsorig : fsjump : epochrng(2)*fsorig);  % NOTE: the 10 makes things nice and even, and hits 0
        ix_raw = onsets(ttt,zzz,p) + (epochrng(1)*fsorig : epochrng(2)*fsorig);	% channel analysis	 (t = -0.5:3.5 s; f = 2 kHz)
        on_samples  = length(epochrng(1)*fsorig:0*fsorig) : length(epochrng(1)*fsorig:0*fsorig) + 2 * fsorig - 1;  % samples corresponding to on stimuli
        %off_samples = cat (2, 1 : length(epochrng(1)*fsorig:0*fsorig)-1, length(epochrng(1)*fsorig:0*fsorig) + 2 * fsorig : 4*fsorig);  % samples corresponding to off-stimuli-period
        off_samples = 1 : length(epochrng(1)*fsorig:0*fsorig)-1;
        temp_bb = bb(sum(bblengths(1:p-1)) + ix_bb);			% bb analysis
		temp_raw = collectdata(sum(bblengths(1:p-1)) + ix_raw);	% raw analysis
        if any(isnan(temp_bb))
          fprintf('BADDATA: ttt=%d, p=%d, ccc=%d, zzz=%d\n',ttt,p,ccc,zzz);
          temp_bb(:) = badval;
        end
        if any(isnan(temp_raw))
		  temp_raw(:) = badval;
        end
        
		%stimco = stimcounter(zzz,ccc,mod2(p,2),a1.stimclassrec(ttt)) + 1;
        bb_data(:,zzz,ccc,mod2(p,2), ...
             stimcounter(zzz,ccc,mod2(p,2),a1.stimclassrec(ttt)) + 1, ...   % a1.stimclassrec tells us the stim number (1-24)
             a1.stimclassrec(ttt)) = temp_bb;
        data(:,zzz,ccc,mod2(p,2), ...
             stimcounter(zzz,ccc,mod2(p,2),a1.stimclassrec(ttt)) + 1, ...   % a1.stimclassrec tells us the stim number (1-24)
            a1.stimclassrec(ttt)) = temp_raw;
        [psdvar,f] = pwelch(temp_raw(on_samples)',hamming(1000),0,2 ^ nextpow2(fsorig),fsorig);
         %figureprep([100 100 900 300],1);plot(f,10*log10(psdvar));
        psd_on(zzz,ccc,mod2(p,2), ...
             stimcounter(zzz,ccc,mod2(p,2),a1.stimclassrec(ttt)) + 1, ...   % calculating psd from raw data
             a1.stimclassrec(ttt),:) = psdvar(1:fupper);
        [psdvar,f] = pwelch(temp_raw(off_samples)',hamming(1000),0,2 ^ nextpow2(fsorig),fsorig);
        psd_off(zzz,ccc,mod2(p,2), ...
             stimcounter(zzz,ccc,mod2(p,2),a1.stimclassrec(ttt)) + 1, ...   % calculating psd from raw data
             a1.stimclassrec(ttt),:) = psdvar(1:fupper);
		stimcounter(zzz,ccc,mod2(p,2),a1.stimclassrec(ttt)) = ...
          stimcounter(zzz,ccc,mod2(p,2),a1.stimclassrec(ttt)) + 1;
      end
    
    end
  end  
end

%% CHANNEL AVG TIME AND FREQUENCY RESPONSE

psd_baseline = squeeze(mean(reshape(psd_off,zzz*numchannels*numtasks*numreps*numstimuli,fupper),1));
bb_data_avg = mean(reshape(bb_data,length(epochtime_bb),numchannels,numtasks*numreps*numstimuli),3)';
data_avg = mean(reshape(data,length(epochtime),numchannels,numtasks*numreps*numstimuli),3)';
psd_avg = squeeze(mean(reshape(psd_on,numchannels,numtasks*numreps*numstimuli,fupper),2));

%% PLOT CHANNEL AVG TIME AND FREQUENCY RESPONSE

clims = [0 50];
figure,imagesc(bb_data_avg,clims);
colorbar;
set(gca,'XTick',0:100:800)
set(gca,'XTickLabel',-0.5:0.5:3.5)

figure,imagesc(data_avg);
colorbar;
set(gca,'XTick',0:1000:8000)
set(gca,'XTickLabel',-0.5:0.5:3.5)

plot(data_avg(75,:))
hold on;
plot(data_avg(76,:))
plot(data_avg(77,:))
plot(data_avg(78,:))
hold off;

%% Plot normalized PSD

ccc = 27;
plot((normalize(10*log10(psd_avg(ccc,:)) - 10*log10(psd_baseline))))

%% Stimuli Timecourse

ccc = 75;

subplot(2,4,1)
plot((squeeze(mean(bb_data(:,1,ccc,1,:,6),5))));
ylim([-50 350]); xlim([0 800]); hold on; 
plot((squeeze(mean(bb_data(:,1,ccc,1,:,7),5))));
plot((squeeze(mean(bb_data(:,1,ccc,1,:,8),5))));
plot((squeeze(mean(bb_data(:,1,ccc,1,:,9),5))));
plot((squeeze(mean(bb_data(:,1,ccc,1,:,4),5))));
hold off;
title('WP_F');

subplot(2,4,5)
plot((squeeze(mean(bb_data(:,1,ccc,2,:,6),5))));
ylim([-50 350]); xlim([0 800]); hold on; 
plot((squeeze(mean(bb_data(:,1,ccc,2,:,7),5))));
plot((squeeze(mean(bb_data(:,1,ccc,2,:,8),5))));
plot((squeeze(mean(bb_data(:,1,ccc,2,:,9),5))));
plot((squeeze(mean(bb_data(:,1,ccc,2,:,4),5))));
hold off;
title('WP_C');

subplot(2,4,2)
plot((squeeze(mean(bb_data(:,1,ccc,1,:,10),5))));
ylim([-50 350]); xlim([0 800]); hold on; 
plot((squeeze(mean(bb_data(:,1,ccc,1,:,11),5))));
plot((squeeze(mean(bb_data(:,1,ccc,1,:,12),5))));
plot((squeeze(mean(bb_data(:,1,ccc,1,:,13),5))));
plot((squeeze(mean(bb_data(:,1,ccc,1,:,5),5))));
hold off;
title('FP_F');

subplot(2,4,6)
plot((squeeze(mean(bb_data(:,1,ccc,2,:,10),5))));
ylim([-50 350]); xlim([0 800]); hold on; 
plot((squeeze(mean(bb_data(:,1,ccc,2,:,11),5))));
plot((squeeze(mean(bb_data(:,1,ccc,2,:,12),5))));
plot((squeeze(mean(bb_data(:,1,ccc,2,:,13),5))));
plot((squeeze(mean(bb_data(:,1,ccc,2,:,5),5))));
hold off;
title('FP_C');

subplot(2,4,3)
plot((squeeze(mean(bb_data(:,1,ccc,1,:,14),5))));
ylim([-50 350]); xlim([0 800]); hold on; 
plot((squeeze(mean(bb_data(:,1,ccc,1,:,15),5))));
plot((squeeze(mean(bb_data(:,1,ccc,1,:,16),5))));
plot((squeeze(mean(bb_data(:,1,ccc,1,:,4),5))));
hold off;
title('WC_F');

subplot(2,4,7)
plot((squeeze(mean(bb_data(:,1,ccc,2,:,14),5))));
ylim([-50 350]); xlim([0 800]); hold on;  
plot((squeeze(mean(bb_data(:,1,ccc,2,:,15),5))));
plot((squeeze(mean(bb_data(:,1,ccc,2,:,16),5))));
plot((squeeze(mean(bb_data(:,1,ccc,2,:,4),5))));
hold off;
title('WC_C');

subplot(2,4,4)
plot((squeeze(mean(bb_data(:,1,ccc,1,:,17),5))));
ylim([-50 350]); xlim([0 800]); hold on; 
plot((squeeze(mean(bb_data(:,1,ccc,1,:,18),5))));
plot((squeeze(mean(bb_data(:,1,ccc,1,:,19),5))));
plot((squeeze(mean(bb_data(:,1,ccc,1,:,5),5))));
hold off;
title('FC_F');

subplot(2,4,8)
plot((squeeze(mean(bb_data(:,1,ccc,2,:,17),5))));
ylim([-50 350]); xlim([0 800]); hold on; 
plot((squeeze(mean(bb_data(:,1,ccc,2,:,18),5))));
plot((squeeze(mean(bb_data(:,1,ccc,2,:,19),5))));
plot((squeeze(mean(bb_data(:,1,ccc,2,:,5),5))));
hold off;
title('FC_C');
