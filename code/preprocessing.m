%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETUP

%path Documents\git_repo\readingECoG\code

addpath(genpath('/Users/zeeshanqadir/Documents/git_repo/knkutils'))

% define
subjects = {'YBA' 'YBD'};
subjectfilenums = {[12 13 14 15] [30 31 32 33]};  % this is for the file names (one file per run per [subject])
channels = 1:128;
numchannels = 128;
photodiode = 129;   % photodiode signal is in #129
epochrng = [-.5 4.5];  % pull out this time range (in seconds)
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
spectratime = (epochrng(1)*fsorig : (epochrng(2)-1)*fsorig)/fsorig; % for spectra to reduce size (-0.5 to 2.5s)
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETUP CONTINUED (ASSIGN LABELS TO CHANNELS)

% define
cfun = @(x,y) cellfun(@(z) [x z],mat2cellstr(y),'UniformOutput',0); % mat2cellstr: return a cell matrix of strings (knkutils)

% initialize
channellabels = {};  % 128 channels for each subject!

% YBA:
channellabels{1} = [ ...
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


% YBD:
channellabels{2} = [ ...
cfun('LF',1:4) ...
cfun('LAT',1:6) ...
cfun('LAST',1:4) ...
cfun('LMST',1:4) ...
cfun('LPST',1:6) ...
cfun('LP',1:6) ...
cfun('LD',1:8) ...
cfun('LTO',1:8) ...
cfun('LMini',1:12) ...
repmat({'nolabel'},[1 6]) ...
cfun('RF',1:4) ...
cfun('RAT',1:4) ...
cfun('RAST',1:6) ...
cfun('RMST',1:4) ...
cfun('RPST',1:6) ...
cfun('RP',1:6) ...
cfun('RD',1:8) ...
cfun('RTO',1:8) ...
cfun('RMini',1:12) ...
repmat({'nolabel'},[1 4]) ...
cfun('EKG',1:2) ...
];

% initialize
anatlabels = {};  % 128 channels for each subject!

% YBA:
anatlabels{1} = [ ...
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

% YBD:
anatlabels{2} = [ ...
repmat({'F'}, [1 4]) ...
repmat({'AT'},[1 6]) ...
repmat({'AT'},[1 4]) ...
repmat({'AT'},[1 4]) ...
repmat({'IT'},[1 6]) ...
repmat({'P'},[1 6]) ...
repmat({'DD'},[1 8]) ...
repmat({'O'},[1 8]) ...
repmat({'M'},[1 12]) ...
repmat({'X'},[1 6]) ...
repmat({'F'},[1 4]) ...
repmat({'AT'},[1 4]) ...
repmat({'AT'},[1 6]) ...
repmat({'AT'},[1 4]) ...
repmat({'IT'},[1 6]) ...
repmat({'P'},[1 6]) ...
repmat({'DD'},[1 8]) ...
repmat({'O'},[1 8]) ...
repmat({'M'},[1 12]) ...
repmat({'X'},[1 4]) ...
repmat({'X'},[1 2]) ...
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

%hacks below - data size > 64GB so break across task + spectra time window
%(-0.5:2.5 s)
% spectra_fix = zeros(length(spectratime),length(subjects),numchannels,numreps,numstimuli,77,'single'); %check how to get 77??
% spectra_cat = zeros(length(spectratime),length(subjects),numchannels,numreps,numstimuli,77,'single');

% psd_on = zeros(length(subjects),numchannels,numtasks,numreps,numstimuli,fupper,'single');
% psd_off = zeros(length(subjects),numchannels,numtasks,numreps,numstimuli,fupper,'single');


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
    figureprep([100 100 1000 300]); hold on;		% figureprep & figurewrite(knkutils)
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
  for ccc= 1:numchannels
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
    
    % compute moving average to reduce data rate to 200Hz (2000/10) in the
    % next step
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
        
        temp_raw(isnan(temp_raw)) = 0;
        [S,f_spectra] = getWaveletSpectrogram(temp_raw, fsorig, [1, fupper]);   % Returns the Morlet (Gabor) wavelet transform (Spectrogram) for a signal - HH
        
        %Uncomment below for spectral analysis (This requies a huge
        %memory!!)
        
%         if mod2(p,2) == 1
%             spectra_fix(:,zzz,ccc,stimco,a1.stimclassrec(ttt),:) = S(:,1:length(spectratime))';
%         elseif mod2(p,2) == 2
%             spectra_cat(:,zzz,ccc,stimco,a1.stimclassrec(ttt),:) = S(:,1:length(spectratime))';
%         end
        
        [psdvar,f] = pwelch(temp_raw(on_samples)',hamming(1000),0,2 ^ nextpow2(fsorig),fsorig);
         %figureprep([100 100 900 300],1);plot(f,10*log10(psdvar));
        psd_on(zzz,ccc,mod2(p,2), stimco, ...
             a1.stimclassrec(ttt),:) = psdvar(1:fupper);
        [psdvar,f] = pwelch(temp_raw(off_samples)',hamming(1000),0,2 ^ nextpow2(fsorig),fsorig);
        psd_off(zzz,ccc,mod2(p,2), stimco, ...
             a1.stimclassrec(ttt),:) = psdvar(1:fupper);
		
        stimcounter(zzz,ccc,mod2(p,2),a1.stimclassrec(ttt)) = ...
           stimcounter(zzz,ccc,mod2(p,2),a1.stimclassrec(ttt)) + 1;
      end
    
    end
  end  
end

%%  Baseline subtraction and normalization

% save('Feb21_exttime_psd.mat','-v7.3');
% load('Feb21_exttime_psd.mat');

% new: mean across time -500 to 0 ms, then across 6 trials, all stim, and then 2 tasks
bb_base = nanmean( nanmean( nanmean( nanmean( bb_data(1:100, :, :, :, :, :), 1), 5), 6), 4);
bbdata_pc = bsxfun(@rdivide, bb_data, bb_base) - 1;

% old
% bb_base = nanmean(nanmean(bb_data(1:100, :, :, :, :, :),5),1); % mean across trials, then mean across time 0:500ms
% bbdata_br = bb_data - nanmean(nanmean(bb_data(1:100, :, :, :, :, :),5),1); % mean across trials, then mean across time 0:500ms
% bbdata_pc = bsxfun(@rdivide,bbdata_br,nanmean(nanmean(bb_base(:, :, :, :, :, setdiff(1:24,[1 3])),6),4));

% take log10 here:
% also take mean baseline across everything (fix + cat) first
% spectra_fix = log(spectra_fix) - log(mean(spectra_fix(1:1001,zzz,ccc,stimco,a1.stimclassrec(ttt),:),1));
% spectra_cat = log(spectra_cat) - log(mean(spectra_cat(1:1001,zzz,ccc,stimco,a1.stimclassrec(ttt),:),1));

%plotting spectra
%load('dg_colormap.mat', 'cm') %KM color scheme;
%figure,uimagesc(spectratime,f_spectra,squeeze(spectra_cat(:,1,75,1,5,:))',[-5 5]); axis xy; colormap(cm); colorbar

%% Stats

% task_number = {[1 2 3 4]   [5  6  7  8]  [9  10 11 12] [13 14 15 16]};
% stimgroups  = {[6 7 8 9]   [10 11 12 13] [14 15 16 4]  [17 18 19 5]};% [20 21 22 10] [2 23 24]};
% stimleg  = {["0", "25", "50", "75"]  ["0", "25", "50", "75"] ["3" "5" "8" "100"] ["4" "6" "10" "100"]};
% stimgrnames = {'Word Phase' 'Face Phase'    'Word Contrast'   'Face Contrast'};%   'Noise Con'   'Other'};


stimgroups  = {4, 5};
stimleg  = {"Word", "Face"};
stimgrnames = {'Word' 'Face'};

p = zeros(length(subjects), numchannels);
t = zeros(length(subjects), numchannels);
r = zeros(length(subjects), numchannels);
randstat = zeros(length(subjects), numchannels);

for zzz = 1:2
    for ccc = 1:numchannels
        x = [];
        for ppp = 1:2
            for stimg = 1:2
                for stimt = 1:numreps
                     x = [x; squeeze(mean(bbdata_pc(121:171, zzz, ccc, ppp, stimt, stimgroups{stimg}), 1))]; % 100 ms - 350 ms
                end
            end
        end
        [~, pval, ~, stats] = ttest(x);
        p (zzz, ccc) = pval;
        t (zzz, ccc) = stats.tstat;
%         d (zzz, ccc) = mean(x)/ sqrt(var(x));
        r (zzz, ccc) = corrcoef(x)^2;
        randstat (zzz, ccc)  = mean(x);
    end
end

[p_fdr, p_th] = ccepPCC_fdr(p,0.05); % fdr -DH

for zzz = 1 : length(subjects)
    figure,scatter (1 : numchannels, squeeze(p_fdr(zzz, :))');
    hold on
    yline(0.05);
end


%% Good Channels(Visually Identified)

%EVCgcc = {horzcat(74:75),horzcat(43:44, 108)};

%LVCgcc = {horzcat(76:78),horzcat(45:46, 109:110)};

%gcc = sort(horzcat((9:11)*2-1, (44:46)*2-1, (74:78)*2-1, (3)*2, (20:22)*2, (86:87)*2, (108:110)*2));

% gcc = {horzcat((9:11), (44:46), (74:78)), horzcat((3), (20:22), (86:87), (108:110))};

gcc = {horzcat((9:11), 15, 41, (44:46), (74:78)), horzcat((20:22), 28, 44, (86:87), 108, 110) };

%% Category Selectivity

d = zeros(length(subjects), numchannels);
y = zeros(2, 2*numreps);

for zzz = 1:2
    for ccc = gcc{zzz}
        for stimg = 1:2
            x = [];
            for ppp = 1:2
                for stimt = 1:numreps
                     x = [x, squeeze(mean(bbdata_pc(121:171, zzz, ccc, ppp, stimt, stimgroups{stimg}), 1))]; % 100 ms - 350 ms
                end
            end
            y(stimg, :) = x;
        end
        d (zzz, ccc) = (mean(y(1,:)) - mean(y(2,:))) / sqrt(0.5 * (var(y(1,:)) + var(y(2,:))));
    end
end

%% Normalized Power

elec_weights = {};

for zzz = 1:2
    y = [];
    for ccc = gcc{zzz}
        x = [];
        for ppp = 1:2
            for stimg = 1:2
                for stimt = 1:numreps
                     x = [x; squeeze(mean(bbdata_pc(121:171, zzz, ccc, ppp, stimt, stimgroups{stimg}), 1))]; % 100 ms - 350 ms
                end
            end
        end
        y = [y, mean(x)];
    end
    elec_weights{zzz} = normalize(y,'range');
end


%% Categorizing electrode anatomically

% Early Visual
elec_EV = {[74, 75], []};

% Late Visual
elec_LV = {[76, 77, 78], [44, 108, 110]};

% Ventral Temporal

elec_VTC = {horzcat((9:11), 15, 41, (44:46)), horzcat((20:22), 28, (86:87))};

% Category-selective VTC
for zzz = 1:2
    x = [];
    y = [];
    z = [];
    for ccc = elec_VTC{zzz}
        if d (zzz, ccc) >= 0.9
            x = [x, ccc];
        elseif d (zzz, ccc) <= -0.9
            y = [y, ccc]; 
        else
            z = [z, ccc];
        end
    end
    elec_VTC_f {zzz} = y;   % face selective 
    elec_VTC_w {zzz} = x;   % word selective
    elec_VTC_ns{zzz} = z;   %  non-selective
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


%% Stimuli Timecourse

stimgroups  = {[6 7 8 9 4]   [10 11 12 13 5] [14 15 16 4] [17 18 19 5]};% [20 21 22 10] [2 23 24]};
stimleg  = {["0", "25", "50", "75", "100"]  ["0", "25", "50", "75", "100"] ["3" "5" "8" "100"] ["4" "6" "10" "100"]};
stimgrnames = {'Word Phase' 'Face Phase'    'Word Contrast'   'Face Contrast'};%   'Noise Con'   'Other'};
stimresp = zeros(length(epochtime_bb),length(gcc{1})+length(gcc{2}),numtasks,length(stimgroups),5);
%stimresp_c = zeros(length(epochtime_bb),length(subjects),20,numtasks,length(stimgroups),5);
stimgain = zeros(length(gcc{1})+length(gcc{2}),numtasks,length(stimgroups),5);
counter = 1;


% One plot per subject per good electrode
for zzz = 1:length(subjects)
   
%     gcc = [];
%     for ccc = 1 : numchannels
%         if p(zzz, ccc) < 0.05
%             gcc = [gcc, ccc];
%         end
%     end
    
    for ccc = gcc{zzz}
        
        %figureprep([100 100 1700 1100]);
        
        for ttt = 1:2
            
            mx = -Inf; mn = Inf;
            
            for stimg = 1:length(stimgroups)
            
                % plot stimtimecourse for fixation(1)/categorization(2) task
                subplot(4,length(stimgroups),stimg+(ttt-1)*length(stimgroups));
                hold on;
                
                for stimc = 1:length(stimgroups{stimg})

                    stimresp(:,counter,ttt,stimg,stimc) = squeeze(mean(bbdata_pc(:,zzz,ccc,ttt,:,stimgroups{stimg}(stimc)),5));
                    plot(smooth(stimresp(:,counter,ttt,stimg,stimc),40),'color',[0.2*stimc 0.5 1-0.2*stimc]);
              
                end
               
                ax = axis;
                mx = max(mx,max(ax(3:4)));
                mn = min(mn,min(ax(3:4)));
                axis([0 1000 ax(3:4)]);
                set(gca,'XTick',0:100:1000); set(gca,'XTickLabel',-0.5:0.5:4.5);
                xlabel('t (s)');
                ylabel('BB response (x fold)');
                if ttt == 1
                    title(sprintf('%s Fixation',stimgrnames{stimg}));
                elseif ttt == 2
                    title(sprintf('%s Categorization',stimgrnames{stimg}));
                end
                legend(stimleg{stimg});
                
            end
            
        end
        
        for ttt = 1:2
            
            for stimg = 1:length(stimgroups)
                
                subplot(4,length(stimgroups),stimg+(ttt-1)*length(stimgroups));
                hold on;
                ax = axis;
                axis([ax(1:2) mn mx]);
                
            end
            
        end
        
        for stimg = 1:length(stimgroups)

            mx = -Inf; mn = Inf;
            % plot stimtimecourse for gain F v/s C
            subplot(4,length(stimgroups),2*length(stimgroups)+stimg);
            hold on;

            stimgain(counter,:,stimg,1) = mean(mean(bbdata_pc(101:300,zzz,ccc,:,:,stimgroups{stimg}(1)),5),1);
            for stimc = 2:length(stimgroups{stimg})
                
                stimgain(counter,:,stimg,stimc) = mean(mean(bbdata_pc(101:300,zzz,ccc,:,:,stimgroups{stimg}(stimc)),5),1);
                plot([stimc-1,stimc],[stimgain(counter,1,stimg,stimc-1), stimgain(counter,1,stimg,stimc)],'r-', 'LineWidth', 2);
                plot([stimc-1,stimc],[stimgain(counter,2,stimg,stimc-1), stimgain(counter,2,stimg,stimc)],'b-', 'LineWidth', 2);
                errorbar(stimc-1,stimgain(counter,1,stimg,stimc-1),2*std(mean(bbdata_pc(101:300,zzz,ccc,1,:,stimgroups{stimg}(stimc-1)),1),0,5)./sqrt(6),'Color','r');  
                errorbar(stimc-1,stimgain(counter,2,stimg,stimc-1),2*std(mean(bbdata_pc(101:300,zzz,ccc,2,:,stimgroups{stimg}(stimc-1)),1),0,5)./sqrt(6),'Color','b'); 
                
            end
            errorbar(stimc,stimgain(counter,1,stimg,stimc),2*std(mean(bbdata_pc(101:300,zzz,ccc,1,:,stimgroups{stimg}(stimc)),1),0,5)./sqrt(6),'Color','r');  
            errorbar(stimc,stimgain(counter,2,stimg,stimc),2*std(mean(bbdata_pc(101:300,zzz,ccc,2,:,stimgroups{stimg}(stimc)),1),0,5)./sqrt(6),'Color','b');
            
            ax = axis;
            mx = max(mx,max(ax(3:4)));
            mn = min(mn,min(ax(3:4)));
            ylabel('BB response (% change)');
            set(gca,'XTick',1:length(stimgroups{stimg})); set(gca,'XTickLabel',stimleg{stimg});
            title(sprintf('%s Gain',stimgrnames{stimg}));
            legend('F','C')
            
        end
        
        for stimg = 1:length(stimgroups)
                
            subplot(4,length(stimgroups),2*length(stimgroups)+stimg);
            hold on;
            ax = axis;
            axis([ax(1)-1 ax(2)+1 mn mx]);
                
        end
        
        for stimg = 1:length(stimgroups)

            mx = -Inf; mn = Inf;

            %plot reaction times
            subplot(4,length(stimgroups),3*length(stimgroups)+stimg);
            hold on;
            
            for stimc = 1:length(stimgroups{stimg})
                bar(stimc,nanmean(reactiontime(zzz,:,stimgroups{stimg}(stimc)),2));
                errorbar(stimc,nanmean(reactiontime(zzz,:,stimgroups{stimg}(stimc)),2),2*nanstd(reactiontime(zzz,:,stimgroups{stimg}(stimc)),0,2)./sqrt(6),'Color','k','HandleVisibility','off');
            end
            ax = axis;
            mx = max(mx,max(ax(3:4)));
            mn = min(mn,min(ax(3:4)));
            ylabel('Response Time (ms)');
            set(gca,'XTick',1:length(stimgroups{stimg})); set(gca,'XTickLabel',stimleg{stimg});
            title(sprintf('%s Reaction Time',stimgrnames{stimg}));
            
        end
        
        for stimg = 1:length(stimgroups)
                
            subplot(4,length(stimgroups),3*length(stimgroups)+stimg);
            hold on;
            ax = axis;
            axis([ax(1:2) mn mx]);
                
        end
        
        %figurewrite(sprintf('Subj%d_ch%03d_%s',zzz,ccc,channellabels{zzz}{ccc}),-1,[],'stimtimecourse');
        %exportgraphics(gcf,[sprintf('stimtimecourse/Subj%d_ch%03d_%s.eps',zzz,ccc,channellabels{zzz}{ccc})]);
        
        counter = counter+1;
        
    end
    
end

%% F v/s C

temp = [];
temp2 = [];

for stimg = 1:length(stimgroups)
    
    figureprep([100 100 1700 1100]);
    mx = -Inf; mn = Inf;
    for stimc = 1:length(stimgroups{stimg})
        
        subplot(3,length(stimgroups{stimg}),stimc); hold on;
        
        nTrials = size(stimresp(:,:,1,stimg,stimc),2);
        yMean = smooth(mean(stimresp(:,:,1,stimg,stimc),2),40);
        ySEM  = smooth(std(stimresp(:,:,1,stimg,stimc),0,2),40)/(6*sqrt(nTrials));
        ts = tinv([0.025, 0.975], nTrials - 1); % n-1 degrees of freedom -HH
        conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
        plot((yMean),'b-', 'LineWidth', 2);
        fill([1:length(epochtime_bb) length(epochtime_bb):-1:1], conf, [0 0 .5],'EdgeColor',[0 0 .5],'FaceAlpha',.5,'HandleVisibility','off');
        
        
        yMean = smooth(mean(stimresp(:,:,2,stimg,stimc),2),40);
        ySEM  = smooth(std(stimresp(:,:,2,stimg,stimc),0,2),40)/(6*sqrt(nTrials));
        ts = tinv([0.025, 0.975], nTrials - 1); % n-1 degrees of freedom
        conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
        plot((yMean),'r-', 'LineWidth', 2);
        fill([1:length(epochtime_bb) length(epochtime_bb):-1:1], conf, [.5 0 0],'EdgeColor',[.5 0 0],'FaceAlpha',.5,'HandleVisibility','off');
        
        ax = axis;
        mx = max(mx,max(ax(3:4)));
        mn = min(mn,min(ax(3:4)));
        xlim([0 400]);
        ylabel('BB response (% change)');
        xlabel('t (s)');
        set(gca,'XTick',0:100:400); set(gca,'XTickLabel',-0.5:0.5:1.5);
        title(stimgrnames{stimg}+" "+stimleg{stimg}(stimc));
        legend('F','C');
        
    end
    
    for stimc = 1:length(stimgroups{stimg})
        
        subplot(3,length(stimgroups{stimg}),stimc); hold on;
        ax = axis;
        axis([ax(1:2) mn mx]);
        %exportgraphics(gca,[sprintf('FvsC/%s%s/',stimgrnames{stimg},stimleg{stimg}(stimc)) 'STC.eps']);
    
    end
        
    mx = -Inf; mn = Inf;
    
    for stimc = 1:length(stimgroups{stimg})
        
        subplot(3,length(stimgroups{stimg}),length(stimgroups{stimg})+stimc); hold on;
        scaling_factor = squeeze(stimresp(:,:,2,stimg,stimc)) - squeeze(stimresp(:,:,1,stimg,stimc));
        nTrials = size(scaling_factor,2);
        yMean = smooth(mean(scaling_factor,2),40);
        %temp = [temp yMean];
        ySEM  = smooth(std(scaling_factor,0,2),40)/(6*sqrt(nTrials));
        ts = tinv([0.025, 0.975], nTrials - 1); % n-1 degrees of freedom -HH
        conf = [(yMean+ts(1)*ySEM)', (yMean(end:-1:1)+ts(2)*ySEM(end:-1:1))'];
        plot(yMean,'k-', 'LineWidth', 2);
        fill([1:length(yMean) length(yMean):-1:1],conf, [.5 0.5 0.5],'EdgeColor',[.5 0.5 0.5],'FaceAlpha',.5,'HandleVisibility','off');
        %[h,p,ci,stats] = ttest2(mean(stimresp_c(:,:,2,stimg,stimc),2),mean(stimresp(:,:,2,stimg,stimc),2));
        for t = 1:length(yMean)     % CoCoSys Lab
            [~, pval(t)] = ttest(squeeze(stimresp(t,:,2,stimg,stimc)), squeeze(stimresp(t,:,1,stimg,stimc)));
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
        xlim([100 300]);
        ylabel('BB response (% change)');
        xlabel('t (s)');
        set(gca,'XTick',100:100:300); set(gca,'XTickLabel',0:0.5:1);
        %title(stimgrnames{stimg}+" "+stimleg{stimg}(stimc));
        title('Scaling Factor');
        
    end
    
    for stimc = 1:length(stimgroups{stimg})
        
        subplot(3,length(stimgroups{stimg}),length(stimgroups{stimg})+stimc); hold on;
        ax = axis;
        axis([ax(1:2) mn mx]);
        %exportgraphics(gca,[sprintf('FvsC/%s%s/',stimgrnames{stimg},stimleg{stimg}(stimc)) 'SF.eps']);
    
    end
    
    subplot(3,length(stimgroups{stimg}),2*length(stimgroups{stimg})+1); hold on;
    
    for stimc = 1:length(stimgroups{stimg})-1
        
        plot([stimc,stimc+1],[mean(stimgain(:,1,stimg,stimc),1),mean(stimgain(:,1,stimg,stimc+1),1)],'r-', 'LineWidth', 2);
        plot([stimc,stimc+1],[mean(stimgain(:,2,stimg,stimc),1),mean(stimgain(:,2,stimg,stimc+1),1)],'b-', 'LineWidth', 2);
        errorbar(stimc,mean(stimgain(:,1,stimg,stimc),1),2*std(stimgain(:,1,stimg,stimc),0,1)./sqrt(20),'Color','r');  
        errorbar(stimc,mean(stimgain(:,2,stimg,stimc),1),2*std(stimgain(:,2,stimg,stimc),0,1)./sqrt(20),'Color','b'); 
        
    end
    errorbar(stimc+1,mean(stimgain(:,1,stimg,stimc+1),1),2*std(stimgain(:,1,stimg,stimc+1),0,1)./sqrt(20),'Color','r');  
    errorbar(stimc+1,mean(stimgain(:,2,stimg,stimc+1),1),2*std(stimgain(:,2,stimg,stimc+1),0,1)./sqrt(20),'Color','b'); 
    ylabel('BB response (% change)');
    set(gca,'XTick',1:length(stimgroups{stimg})); set(gca,'XTickLabel',stimleg{stimg});
    ax = axis;
    axis([ax(1)-1 ax(2)+1 ax(3:4)]);
    %title(sprintf('%s Gain',stimgrnames{stimg}));
    title('Gain');
    legend('F','C')
    %exportgraphics(gca,[sprintf('FvsC/%s%s/',stimgrnames{stimg},stimleg{stimg}(stimc)) 'Gain.eps']);
    
%     subplot(3,length(stimgroups),2*length(stimgroups{stimg})+2); hold on;
    
    for stimc = 1:length(stimgroups{stimg})
    
        bar(stimc,mean(nanmean(reactiontime(:,:,stimgroups{stimg}(stimc)),2),1));
        temp2 = [temp2 mean(nanmean(reactiontime(:,:,stimgroups{stimg}(stimc)),2),1)];
        errorbar(stimc,mean(nanmean(reactiontime(:,:,stimgroups{stimg}(stimc)),2),1),2*std(nanmean(reactiontime(:,:,stimgroups{stimg}(stimc)),2),0,1),'Color','k','HandleVisibility','off');
      
    end
%     hold off;
%     ylabel('Response Time (ms)');
%     set(gca,'XTick',1:length(stimgroups{stimg})); set(gca,'XTickLabel',stimleg{stimg});
%     title(sprintf('%s Reaction Time',stimgrnames{stimg}));
%     legend(stimleg{stimg});
    %exportgraphics(gcf,sprintf('%s.eps',stimgrnames{stimg}))
    figurewrite(sprintf('%s',stimgrnames{stimg}),[],[],'FvsC');
end

%%

plot(temp2',mean(temp,1)');
scatter(temp2',mean(temp,1)')


%% SF for VWFA and FFA  (overall)

vwfaelec = {[],[]};
ffaelec = {[],[]};
targetstim = {[ 8 9 4 ] [ 12 13 5]};
otherstim = {[ 12 13 5  2 23 24] [8 9 4 2 23 24]};
for zzz = 1:length(subjects)
    for ccc = gcc{zzz}
        for stimc = 1:length(targetstim)
            targetstimresp = squeeze(mean(mean(mean(bbdata_br(101:300,zzz,ccc,2,:,targetstim{stimc}),5),6),1));
            otherstimresp = squeeze(mean(mean(mean(bbdata_br(101:300,zzz,ccc,2,:,otherstim{stimc}),5),6),1));
            if stimc == 1 && targetstimresp > otherstimresp
                vwfaelec{zzz} = horzcat(vwfaelec{zzz}, ccc);
            end
            if stimc == 2 && targetstimresp > otherstimresp
                ffaelec{zzz} = horzcat(ffaelec{zzz}, ccc);
            end 
        end
    end
end

color = { 'r', 'g', 'b', 'y', 'c'};
figureprep([100 100 1700 1100]);
counter =1;
for stimg = 1:2:length(stimgroups)
    subplot(1,2,counter);
    for stimc = 1:length(stimgroups{stimg})
        for zzz = 1: length(subjects)
            vwfarespf{zzz} = squeeze((mean(stimgain(zzz,vwfaelec{zzz},1,stimg,stimc),2)));
            vwfarespc{zzz} = squeeze((mean(stimgain(zzz,vwfaelec{zzz},2,stimg,stimc),2)));
            ffarespf{zzz}  = squeeze((mean(stimgain(zzz, ffaelec{zzz},1,stimg,stimc),2)));
            ffarespc{zzz}  = squeeze((mean(stimgain(zzz, ffaelec{zzz},2,stimg,stimc),2)));
        end
        quiver(mean(cellfun(@mean,vwfarespf)),mean(cellfun(@mean,ffarespf)),mean(cellfun(@mean,vwfarespc))-mean(cellfun(@mean,vwfarespf)),mean(cellfun(@mean,ffarespc))-mean(cellfun(@mean,ffarespf)),0,color{stimc});
        
        hold on;
        for zzz = 1: length(subjects)
            vwfarespf{zzz} = squeeze((mean(stimgain(zzz,vwfaelec{zzz},1,stimg+1,stimc),2)));
            vwfarespc{zzz} = squeeze((mean(stimgain(zzz,vwfaelec{zzz},2,stimg+1,stimc),2)));
            ffarespf{zzz}  = squeeze((mean(stimgain(zzz, ffaelec{zzz},1,stimg+1,stimc),2)));
            ffarespc{zzz}  = squeeze((mean(stimgain(zzz, ffaelec{zzz},2,stimg+1,stimc),2)));
        end
        quiver(mean(cellfun(@mean,vwfarespf)),mean(cellfun(@mean,ffarespf)),mean(cellfun(@mean,vwfarespc))-mean(cellfun(@mean,vwfarespf)),mean(cellfun(@mean,ffarespc))-mean(cellfun(@mean,ffarespf)),0,color{stimc});
        
        hold on;
    end
    hold off;
    legend(stimleg{stimg});
    title(sprintf('SF_%s',stimgrnames{stimg}));
    counter = counter+1;
end
figurewrite('plot',[],[],'SF');

%% SF for VWFA and FFA  (cumulative)

color = { 'r', 'g', 'b', 'm', 'c'};



for inc = 1:length(10:10:200)
    counter = 1;
    for stimg = 1:2:length(stimgroups)
        subplot(1,2,counter);
        for stimc = 1:length(stimgroups{stimg})
            for zzz = 1: length(subjects)
                vwfarespf{zzz} = squeeze(mean(mean(mean(bbdata_pc(101:(101+inc*10),zzz,vwfaelec{zzz},1,:,stimgroups{stimg}(stimc)),5),1),2));
                vwfarespc{zzz} = squeeze((mean(mean(mean(bbdata_pc(101:(101+inc*10),zzz,vwfaelec{zzz},2,:,stimgroups{stimg}(stimc)),5),1),2)));
                ffarespf{zzz}  = squeeze((mean(mean(mean(bbdata_pc(101:(101+inc*10),zzz,ffaelec{zzz},1,:,stimgroups{stimg}(stimc)),5),1),2)));
                ffarespc{zzz}  = squeeze((mean(mean(mean(bbdata_pc(101:(101+inc*10),zzz,ffaelec{zzz},2,:,stimgroups{stimg}(stimc)),5),1),2)));
            end
            quiver(mean(cellfun(@mean,vwfarespf)),mean(cellfun(@mean,ffarespf)),mean(cellfun(@mean,vwfarespc))-mean(cellfun(@mean,vwfarespf)),mean(cellfun(@mean,ffarespc))-mean(cellfun(@mean,ffarespf)),0,color{stimc});
        
            hold on;
            for zzz = 1: length(subjects)
                vwfarespf{zzz} = squeeze((mean(mean(mean(bbdata_pc(101:(101+inc*10),zzz,vwfaelec{zzz},1,:,stimgroups{stimg+1}(stimc)),5),1),2)));
                vwfarespc{zzz} = squeeze((mean(mean(mean(bbdata_pc(101:(101+inc*10),zzz,vwfaelec{zzz},2,:,stimgroups{stimg+1}(stimc)),5),1),2)));
                ffarespf{zzz}  = squeeze((mean(mean(mean(bbdata_pc(101:(101+inc*10),zzz,ffaelec{zzz},1,:,stimgroups{stimg+1}(stimc)),5),1),2)));
                ffarespc{zzz}  = squeeze((mean(mean(mean(bbdata_pc(101:(101+inc*10),zzz,ffaelec{zzz},2,:,stimgroups{stimg+1}(stimc)),5),1),2)));
            end
            quiver(mean(cellfun(@mean,vwfarespf)),mean(cellfun(@mean,ffarespf)),mean(cellfun(@mean,vwfarespc))-mean(cellfun(@mean,vwfarespf)),mean(cellfun(@mean,ffarespc))-mean(cellfun(@mean,ffarespf)),0,color{stimc});
        
            hold on;
        end
        hold off;
        legend(stimleg{stimg});
        title(sprintf('SF%s',stimgrnames{stimg}));
        counter = counter+1;
    end
    figurewrite(sprintf('%04d',inc*50),[],[],'cSF');
end



%% SF for VWFA and FFA  (progressive)


color = { 'r', 'g', 'b', 'm', 'c'};
ttl = {'Phase','','Contrast'};

x01 = 0;
x02 = 0;
x03 = 0;
x04 = 0;
y01 = 0;
y02 = 0;
y03 = 0;
y04 = 0;


for stimg = 1:2:length(stimgroups)
    %subplot(1,2,counter);
    counter = 1;
    for stimc = 1:length(stimgroups{stimg})
        for inc = 1:length(10:10:200)
            for zzz = 1: length(subjects)
                vwfarespf{zzz} = squeeze(mean(mean(mean(bbdata_pc(101:(101+inc*10),zzz,vwfaelec{zzz},1,:,stimgroups{stimg}(stimc)),5),1),2));
                vwfarespc{zzz} = squeeze((mean(mean(mean(bbdata_pc(101:(101+inc*10),zzz,vwfaelec{zzz},2,:,stimgroups{stimg}(stimc)),5),1),2)));
                ffarespf{zzz}  = squeeze((mean(mean(mean(bbdata_pc(101:(101+inc*10),zzz,ffaelec{zzz},1,:,stimgroups{stimg}(stimc)),5),1),2)));
                ffarespc{zzz}  = squeeze((mean(mean(mean(bbdata_pc(101:(101+inc*10),zzz,ffaelec{zzz},2,:,stimgroups{stimg}(stimc)),5),1),2)));
            end
            
            x1 = mean(cellfun(@mean,vwfarespf));
            y1 = mean(cellfun(@mean,ffarespf));
            x2 = mean(cellfun(@mean,vwfarespc));
            y2 = mean(cellfun(@mean,ffarespc));
            quiver(x01,y01,x1-x01,y1-y01,'m');
            hold on;
            quiver(x02,y02,x2-x02,y2-y02,'r');
            hold on;
%             quiver(x1,y1,x2-x1,y2-y1,'k');
%             hold on;
            x01 = x1;
            x02 = x2;
            y01 = y1;
            y02 = y2;
            for zzz = 1: length(subjects)
                vwfarespf{zzz} = squeeze((mean(mean(mean(bbdata_pc(101:(101+inc*10),zzz,vwfaelec{zzz},1,:,stimgroups{stimg+1}(stimc)),5),1),2)));
                vwfarespc{zzz} = squeeze((mean(mean(mean(bbdata_pc(101:(101+inc*10),zzz,vwfaelec{zzz},2,:,stimgroups{stimg+1}(stimc)),5),1),2)));
                ffarespf{zzz}  = squeeze((mean(mean(mean(bbdata_pc(101:(101+inc*10),zzz,ffaelec{zzz},1,:,stimgroups{stimg+1}(stimc)),5),1),2)));
                ffarespc{zzz}  = squeeze((mean(mean(mean(bbdata_pc(101:(101+inc*10),zzz,ffaelec{zzz},2,:,stimgroups{stimg+1}(stimc)),5),1),2)));
            end
            x3 = mean(cellfun(@mean,vwfarespf));
            y3 = mean(cellfun(@mean,ffarespf));
            x4 = mean(cellfun(@mean,vwfarespc));
            y4 = mean(cellfun(@mean,ffarespc));
            quiver(x03,y03,x3-x03,y3-y03,'c');
            hold on;
            quiver(x04,y04,x4-x04,y4-y04,'b');
            hold on;
%             quiver(x3,y3,x4-x3,y4-y3,'g');
%             hold on;
            x03 = x3;
            x04 = x4;
            y03 = y3;
            y04 = y4;
        end
        xlabel('VWFA');
        ylabel('FFA');
        legend(["Word Fix","Word Cat","Face Fix","Face Cat"]);
        title(sprintf('SF-%s-%s',ttl{stimg},stimleg{stimg}(counter)));     
        hold off;
        exportgraphics(gcf,[sprintf('pSF/SF-%s-%s.eps',ttl{stimg},stimleg{stimg}(counter))]);
        figurewrite(sprintf('SF-%s-%s',ttl{stimg},stimleg{stimg}(counter)),[],[],'pSF');
        counter = counter+1;
    end
    
end


%% Plotting Spectograms

tasklabels = {'Fix' 'Cat'};
stimgroups  = {[6 7 8 9]   [10 11 12 13] [14 15 16 4] [17 18 19 5]};% [20 21 22 10] [2 23 24]};
stimleg  = {["0", "25", "50", "75"]  ["0", "25", "50", "75"] ["3" "5" "8" "100"] ["4" "6" "10" "100"]};
stimgrnames = {'Word Phase-coherence' 'Face Phase-coherence'    'Word Contrast'   'Face Contrast'};%   'Noise Con'   'Other'};

spect_base = zeros(length(subjects),numchannels,77);

for zzz = 2%1:2
    for ccc = 108%1:numchannels
        S_temp = [];
        for ppp = 1:2
            for stimg = 1:length(stimgroups)
                for stimc = 1:length(stimgroups{stimg})
                    temp = mean(data(:, zzz, ccc, ppp, :, stimgroups{stimg}(stimc)), 5);
                    temp(isnan(temp)) = 0;
                    [S,f_spectra] = getWaveletSpectrogram(temp, fsorig, [1, fupper]);
                    S_temp = [S_temp; S'];
                end
            end
        end
        spect_base (zzz,ccc,:) = mean(S_temp, 1);
    end
end



for zzz = 2%1:2
    for ccc = 108%1:numchannels
        for ppp = 2%1:2
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
