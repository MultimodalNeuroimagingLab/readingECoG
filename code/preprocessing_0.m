% Script to perfom preprocessing on ECoG data collected on 2 subjects
% 
% Knedrick Kay, Zeeshan Qadir, Dora Hermes 2023

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETUP

dataPath = setLocalDataPath(1);

% define
subjects = {'YBA' 'YBD'};
subjectfilenums = {[12 13 14 15] [30 31 32 33]};  % this is for the file names (one file per run per [subject])
channels        = 1:128;
numchannels     = length(channels);
photodiode      = 129;              % photodiode signal is in #129
epochrng        = [-.5 4.5];        % pull out this time range (in seconds)
onsetix         = {{1:66 1:66 1:66 1:66} {1:66 1:66 1:66 1:66}};  % indices to pull the black circle (based on the no of stimulus trials per run; photodiode detects this circle)
fsorig          = 2000;             % sampling rate (original)
fsjump          = 10;               % moving average length
fs              = fsorig/fsjump;    % sampling rate (down-sampled to)
fupper          = 200;              % upperbound for frequency analysis
numtasks        = 2;                % alternate the tasks (Fixation/Categorization)
numreps         = 6;                % 6 total trials for each numtasks (FC#1: 3 trials, FC#2: 3 trials)
numstimuli      = 24;               % remember that we don't use #1 and #3, so entries for this will be blank 
numtrials       = 66;               % in one run, there are these many stimulus trials (numstimuli (=22) X 3 trials per numstimuli)
numruns         = 4;                % total number of runs per subject (FCFC)
tasklabels      = {'Fixation' 'Categorization'};
 
% Epoch times
epochtime    = (epochrng(1)*fsorig : epochrng(2)*fsorig)/fsorig; % time in seconds for each data point of an epoch
epochtime_bb = (epochrng(1)*fs : epochrng(2)*fs)/fs;             % BB epochtime
spectratime  = find(epochtime >= -0.5 & epochtime < 2.5);        % for spectra to reduce size (-0.5 to 2.5s)
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



%% Finding epochs

onsets = zeros( numtrials, length(subjects), numruns);  % indices of data points at which trials occur
onsetpoints = {};


for zzz=1:length(subjects)

  % get behavioral files
  files = matchfiles( fullfile( dataPath.input, subjects{zzz}, 'LogFiles', '*.mat'))   %matchfiles: to match the filename (knkutils)
  assert( length(files) == length( subjectfilenums{zzz}));
  
  temponsetpoints = []; 
  optemp = 0;
  % process each run (photodiode)
  for p=1:length(subjectfilenums{zzz})

    % load photodiode information
    chantoload = photodiode;     % the special 129
    file0 = fullfile( dataPath.input, subjects{zzz}, ...
                      sprintf('%sDatafile%03d_ch%d.mat', subjects{zzz},subjectfilenums{zzz}(p),chantoload));

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
    figurewrite(sprintf('photodiode_subj%d_file%d',zzz,p),[],[], fullfile( dataPath.output, 'photodiode'));

    % record
    onsets(:,zzz,p) = theseOnsets;
    temponsetpoints = [temponsetpoints optemp+theseOnsets]; % cumulative epoch points
    optemp = optemp + length(pd.analogTraces);	% cumulative length of runs 
  end
  onsetpoints{zzz} = single(temponsetpoints);
 
end

clear zzz files temponsetpoints optemp p chantoload file0 pd fs0 numsam theseOnsets okixs tempdata locoftransition

%%  Channel Analysis

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PREP DATA

% do it

data       = zeros(length(epochtime),length(subjects),numchannels,numtasks,numreps,numstimuli,'single');
bb_data    = zeros(length(epochtime_bb),length(subjects),numchannels,numtasks,numreps,numstimuli,'single');
alpha_data = zeros(length(epochtime_bb),length(subjects),numchannels,numtasks,numreps,numstimuli,'single');

% hacks below - data size > 64GB so break across task + spectra time window (-0.5:2.5 s)
% spectra_fix = zeros(length(spectratime),length(subjects),numchannels,numreps,numstimuli,77,'single'); %check how to get 77??
% spectra_cat = zeros(length(spectratime),length(subjects),numchannels,numreps,numstimuli,77,'single');

psd_on  = zeros(length(subjects),numchannels,numtasks,numreps,numstimuli,fupper,'single');
psd_off = zeros(length(subjects),numchannels,numtasks,numreps,numstimuli,fupper,'single');

stimcounter = zeros(length(subjects),numchannels,numtasks,numstimuli);  % number of trials encountered so far

recdata = {};   % continuous (un-epoched) raw data for each electorde
bbtemp = {};    % continuous (un-epoched) bb time series for each electode

%%%%%%%%

for zzz=1:length(subjects)

  % Ignoring non-recording/non-ECoG/mini channels
  if zzz == 1
      good_channels = [1:26, 29:78];
  elseif zzz ==2
      good_channels = [1:46, 65:110];
  else
      assert(1==0);
  end
    
  % get behavioral files
  files = matchfiles( fullfile( dataPath.input, subjects{zzz}, 'LogFiles', '*.mat'))   %matchfiles: to match the filename (knkutils)
  assert(length(files)==length(subjectfilenums{zzz}));
  
  signal_raw = [];
  bblengths = zeros( numchannels, 4);

  % process each channel
  for ccc= 1:numchannels
    
    % init
    collectdata = [];
    
    % process each run (actual data)
    for p=1:length(subjectfilenums{zzz})

      % load data
      chantoload = channels(ccc);  % the usual 1-128
      file0 = fullfile( dataPath.input, subjects{zzz}, ...
                        sprintf('%sDatafile%03d_ch%d.mat', subjects{zzz},subjectfilenums{zzz}(p),chantoload));
      if ~exist(file0,'file')
        continue;
      end
      pd = load(file0);
      fs0 = pd.analogInfos.SampleRate; assert(fs0==fsorig);
      
      % record data
      collectdata = [collectdata pd.analogTraces];
      bblengths(ccc, p) = length(pd.analogTraces);

    end
    
    % TOTAL HACK SPIKE (detect later!)
    if zzz==1   %subject YBA
      badrng = 1.072 * 10^6 : 1.096 * 10^6; % data seems to be corrupted in this time range for all channels
      %bb(badrng) = NaN;  %%median(bb(setdiff(1:length(bb),badrng)));
      collectdata(badrng) = NaN;
      %collectdata(badrng) = nanmedian(collectdata);
    end

    % Some channels didn't record for entire duration, hence ot fix
    if ~isempty(signal_raw) && (length(collectdata) ~= size(signal_raw,1))
        collectdata( end : size(signal_raw,1)) = nan;
    end

    signal_raw = [signal_raw collectdata'];     % time x channels

  end
    
  temp_signalraw = signal_raw;
  temp_signalraw( isnan( temp_signalraw)) = 0;

  % notch filter for line noise
  signal_notch = ieeg_notch( temp_signalraw, fsorig, 60, [], true);
  
  clear temp_signalraw
    
  % reref - good channels only

  data_out = zeros( size( signal_notch'));

  % CAR in 64-length block
  x = ieeg_carRegress( signal_notch', good_channels( good_channels<65));
  data_out(1:64,:) = x(1:64,:);
  x = ieeg_carRegress( signal_notch', good_channels( good_channels>64 & good_channels<129)); 
  data_out(65:128,:) = x(65:128,:);

  clear x

  data_out = data_out';

  %%% BB analysis

  for ccc= 1:numchannels

    collectdata = data_out( :, ccc);

    recdata{zzz,ccc} = single(collectdata);    % raw data amplitude

    % if no data exists, get out early
    if isempty(collectdata)
      continue;
    end

    sprintf('Subject%d_%03d',zzz,ccc)

	% NEW way:
    bands = [70  90;  % 20 hz bins, avoiding 60 and 180 (max 210)
             90  110;
             110 130;
             130 150;
             150 170];   %% HACK OUT
%             190 210];
    alpha = [8 12];

    %check the bb influence on nb
    
    bb = ecog_extractBroadband(collectdata,fsorig,[],bands);  % NOTE: WE DO WHOLE EXPERIMENT AT ONCE TO ENSURE EQUAL SCALING ISSUE...
                                                              % ecog_extractBroadband: mean after hilbert; power; geomean (ECoG utilities by JW_NYU)
    alphab = ecog_extractBroadband(collectdata,fsorig,[],alpha);
    
    %bb_bp = ieeg_butterpass(collectdata, [70 170], fsorig);   % bandpass; amplitude (mnl_ieegBasics by DH)

	% 	RAW CASE: bb = pd.analogTraces;

    
    bbtemp{zzz,ccc} = single(bb);   % save broadband so we can inspect it!

    
    %assert(all(isfinite(bb)));
    
    % compute moving average to reduce data rate to 200Hz (2000/10) in the
    % next step
    bb = movmean(bb,  fsjump,1);  % slight slop due to evenness of the number.  we can plop this in the photodiode transition issue.
    alphab = movmean(alphab,  fsjump,1);

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
        
        temp_bb = bb(sum(bblengths(ccc, 1:p-1)) + ix_bb);			% bb analysis
        temp_alpha = alphab(sum(bblengths(ccc, 1:p-1)) + ix_bb);     % alpha band analysis

		temp_raw = collectdata(sum(bblengths(ccc, 1:p-1)) + ix_raw);	% raw analysis
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
        alpha_data(:,zzz,ccc,mod2(p,2), stimco, ...
             a1.stimclassrec(ttt)) = temp_alpha;
         
        data(:,zzz,ccc,mod2(p,2), stimco, ...
             a1.stimclassrec(ttt)) = temp_raw;
        
        temp_raw(isnan(temp_raw)) = 0;
        %[spectra,f_spectra] = getWaveletSpectrogram(temp_raw, fsorig, [1, fupper]);   % Returns the Morlet (Gabor) wavelet transform (Spectrogram) for a signal - HH
        
        %Uncomment below for spectral analysis (This requies a huge
        %memory!!)
        
%         if mod2(p,2) == 1
%             spectra_fix(:,zzz,ccc,stimco,a1.stimclassrec(ttt),:) = spectra(:,1:length(spectratime))';
%         elseif mod2(p,2) == 2
%             spectra_cat(:,zzz,ccc,stimco,a1.stimclassrec(ttt),:) = spectra(:,1:length(spectratime))';
%         end
        
        [psdvar,f_psd] = pwelch(temp_raw(on_samples)',hamming(1000),0,2 ^ nextpow2(fsorig),fsorig);
         %figureprep([100 100 900 300],1);plot(f,10*log10(psdvar));
        psd_on(zzz,ccc,mod2(p,2), stimco, ...
             a1.stimclassrec(ttt),:) = psdvar(1:fupper);
        [psdvar,f_psd] = pwelch(temp_raw(off_samples)',hamming(1000),0,2 ^ nextpow2(fsorig),fsorig);
        psd_off(zzz,ccc,mod2(p,2), stimco, ...
             a1.stimclassrec(ttt),:) = psdvar(1:fupper);
        f_psd = f_psd(1:fupper);
		 
        stimcounter(zzz,ccc,mod2(p,2),a1.stimclassrec(ttt)) = ...
           stimcounter(zzz,ccc,mod2(p,2),a1.stimclassrec(ttt)) + 1;
      end
    
    end
  end  
end

clear zzz files ccc collectdata bblengths p chantoload file0 pd fs0 bb alphab a1 ttt temp* psdvar stimco ix* signal_raw signal_notch

save(fullfile( dataPath.output, 'preprocessing_0.mat'), '-v7.3');