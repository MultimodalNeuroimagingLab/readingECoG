%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETUP

% define
subjects = {'YBA'};
subjectfilenums = {[12 13 14 15]};  % this is for the file names (one file per run?)
channels = 1:128;
numchannels = 128;
photodiode = 129;   % photodiode signal is in #129
epochrng = [-.5 3.5];  % pull out this time range (in seconds)
onsetix = {{1:66 1:66 1:66 1:66}};  % indices to pull the black circle (based on the no of stimulus trials per run; photodiode detects this circle)
fsorig = 2000;    % sampling rate (original)
fsjump = 10;      % moving average length
fs = 200;         % sampling rate (down-sampled to)
numtasks = 2;     % alternate the tasks
numreps = 6;      % 6 total trials for each numtasks (FC_1 - 3 trials, FC_2 = 3 trials)
numstimuli = 24;  % remember that we don't use #1 and #3, so entries for this will be blank
numtrials = 66;   % in one run, there are these many stimulus trials (numstimuli X 3/numstimuli)
numruns = 4;      % total number of runs per subject (FCFC)
tasklabels = {'Fixation' 'Categorization'};
 
% calc
epochtime = (epochrng(1)*fs : epochrng(2)*fs)/fs;  % time in seconds for each data point of an epoch (5 ms interval)

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
psd = zeros(length(subjects),numchannels,numtasks,numreps,numstimuli,101,'single');
stimcounter = zeros(length(subjects),numchannels,numtasks,numstimuli);  % number of trials encountered so far
recdata = {};
bbtemp = {};
for zzz=1:length(subjects)

  % get behavioral files
  files = matchfiles(sprintf('%s/LogFiles/*.mat',subjects{zzz}))   %matchfiles: to match the filename (knkutils)
  assert(length(files)==length(subjectfilenums{zzz}));

  % process each run (photodiode)
  for p=1:length(subjectfilenums{zzz})

    % load photodiode information
    chantoload = photodiode;     % the special 129
    file0 = sprintf('%s/%sDatafile%03d_ch%d.mat', ...
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
    figureprep([100 100 1000 300]); hold on;
    plot(pd.analogTraces);
    straightline(theseOnsets,'v','m-');
    figurewrite(sprintf('photodiode_subj%d_file%d',zzz,p),[],[],'~/inout/photodiode');

    % record
    onsets(:,zzz,p) = theseOnsets;
  
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
      file0 = sprintf('%s/%sDatafile%03d_ch%d.mat', ...
                      subjects{zzz},subjects{zzz},subjectfilenums{zzz}(p),chantoload);
      if ~exist(file0,'file')
        continue;
      end
      pd = load(file0);
      fs0 = pd.analogInfos.SampleRate; assert(fs0==fsorig);

% NO LONGER VALID
%       % calc
%       assert(length(pd.analogTraces)==numsam);  % check that same as photodiode length
      
      % record data
      collectdata = [collectdata pd.analogTraces];
      bblengths(p) = length(pd.analogTraces);

    end
    
    % if no data exists, get out early
    if isempty(collectdata)
      continue;
    end
    
    % now do broadband:

% OLD way (version 1):
%        bb = g(pd.analogTraces);   % <===== NOTE THIS!!! 

% NEW way:
    bands = [70 90;  % 20 hz bins, avoiding 60 and 180 (max 210)
             90 110;
             110 130;
             130 150;
             150 170];   %% HACK OUT
%             190 210];
    bb = ecog_extractBroadband(collectdata',fsorig,[],bands);  % NOTE: WE DO WHOLE EXPERIMENT AT ONCE TO ENSURE EQUAL SCALING ISSUE...
                                                               % ecog_extractBroadband (ECoG utilities by JW_NYU)
                                                        
% RAW CASE:
%        bb = pd.analogTraces;

    % save broadband so we can inspect it!
    bbtemp{zzz,ccc} = single(bb);
    recdata{zzz,ccc} = single(collectdata');
    

    % TOTAL HACK SPIKE (detect later!)
    assert(all(isfinite(bb)));
    if zzz==1   %subject YBA
      badrng = 1.072 * 10^6 : 1.09 * 10^6;
      bb(badrng) = NaN;  %%median(bb(setdiff(1:length(bb),badrng)));
      badval = nanmedian(bb);  % **ROBUST
    end
    
    % compute moving average to reduce data rate to 200Hz (2000/10)
    bb = movmean(bb,  fsjump,1);  % slight slop due to evenness of the number.  we can plop this in the photodiode transition issue.
    collectdata = movmean(collectdata',  fsjump,1);
    
    % process each run
    for p=1:length(subjectfilenums{zzz})

      % load behavioral file
      a1 = load(files{p});

      % extract epochs (trials)
      for ttt=1:size(onsets,1)
        ix = onsets(ttt,zzz,p) + (epochrng(1)*fsorig : fsjump : epochrng(2)*fsorig);  % NOTE: the 10 makes things nice and even, and hits 0
        temp = bb(sum(bblengths(1:p-1)) + ix);
        temp1 = collectdata(sum(bblengths(1:p-1)) + ix);
        if any(isnan(temp))
          fprintf('BADDATA: ttt=%d, p=%d, ccc=%d, zzz=%d\n',ttt,p,ccc,zzz);
          temp(:) = badval;
        end
        data(:,zzz,ccc,mod2(p,2), ...
             stimcounter(zzz,ccc,mod2(p,2),a1.stimclassrec(ttt)) + 1, ...   % a1.stimclassrec tells us the stim number (1-24)
             a1.stimclassrec(ttt)) = temp;
        [psdvar,f] = pwelch(temp1,400,0,fs,fs);
        %stimco = stimcounter(zzz,ccc,mod2(p,2),a1.stimclassrec(ttt)) + 1;
        psd(zzz,ccc,mod2(p,2), ...
             stimcounter(zzz,ccc,mod2(p,2),a1.stimclassrec(ttt)) + 1, ...   % calculating psd from raw data
             a1.stimclassrec(ttt),:) = psdvar';
        %figureprep([100 100 900 300],1);plot(f,10*log10(psdvar));
        %figurewrite(sprintf('~/psd/task%d_%03d%02d%d', ...
        %    mod2(p,2),ccc,a1.stimclassrec(ttt),stimco));
        stimcounter(zzz,ccc,mod2(p,2),a1.stimclassrec(ttt)) = ...
          stimcounter(zzz,ccc,mod2(p,2),a1.stimclassrec(ttt)) + 1;
      end

    end
  end
end

%% PSD AVG ACROSS TRIALS FOR A CHANNEL
psdmean_trials= zeros(length(subjects),numchannels,numtasks,numstimuli,length(psdvar),'single');

for lll =1:numstimuli
    psdmean_trials(:,:,:,lll,:) = (psd(:,:,:,1,lll,:)+psd(:,:,:,2,lll,:)+psd(:,:,:,3,lll,:))/3;    
end


tempsum = 0;
for p =1:numtasks
    for lll =1:numstimuli
        tempsum = tempsum + squeeze(psdmean_trials(:,:,p,lll,:));
    end
    for lll =1:numstimuli
        psdmean_trials(:,:,p,lll,:) = log(squeeze(psdmean_trials(:,:,p,lll,:))./(tempsum/numstimuli));
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


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% QUICK LOOK AT SOME BB TIMECOURSES TO CHECK SANITY

figureprep([100 100 900 300],1); plot(log(smooth(bbtemp{1,46},2000)));  %subject YBA
figurewrite('bb_subjYBA_ch46');
figureprep([100 100 900 300],1); plot(log(smooth(bbtemp{1,78},2000)));  %subject YBA
figurewrite('bb_subjYBA_ch78');
%figureprep([100 100 900 300],1); plot(log(smooth(bbtemp{3,86},2000)));
%figureprep([100 100 900 300],1); plot(log(smooth(bbtemp{3,108},2000)));

figure;plot(mean(catcell(2,bbtemp(1,:)),2));
figurewrite('mean');
%figure;plot(mean(catcell(2,bbtemp(3,:)),2));

% OLD:
%         % visualize for sanity
%         figureprep([100 100 1000 300]); hold on;
%         plot(bb);
%         figurewrite(sprintf('ch%d',chantoload),[],[],sprintf('~/inout/bb_subj%d/file%d',zzz,p));

