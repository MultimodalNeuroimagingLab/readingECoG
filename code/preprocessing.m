%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SETUP

%path D:\Zeeshan\reading-ecog\readingECoG

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
%epochtime = (epochrng(1)*fs : epochrng(2)*fs)/fs;
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
psd = zeros(length(subjects),numchannels,numtasks,numreps,numstimuli,fupper,'single');
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
    % figureprep([100 100 1000 300]); hold on;		% figureprep & figurewrite(knkutils)
    figure('Position',[100 100 1000 300]); hold on;
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
        ix = onsets(ttt,zzz,p) + (epochrng(1)*fsorig : fsjump : epochrng(2)*fsorig);  % NOTE: the 10 makes things nice and even, and hits 0
        ix1 = onsets(ttt,zzz,p) + (0*fsorig : 2*fsorig);		% frequency analysis (t =   0 : 2  s; f = 2 kHz)
        ix2 = onsets(ttt,zzz,p) + (-0.5*fsorig : 3.5*fsorig);	% channel analysis	 (t = -0.5:3.5 s; f = 2 kHz)
        temp = bb(sum(bblengths(1:p-1)) + ix);				% bb analysis
        temp1 = collectdata(sum(bblengths(1:p-1)) + ix1);	% frequency analysis
		temp2 = collectdata(sum(bblengths(1:p-1)) + ix2);	% channel analysis
        if any(isnan(temp))
          fprintf('BADDATA: ttt=%d, p=%d, ccc=%d, zzz=%d\n',ttt,p,ccc,zzz);
          temp(:) = badval;
        end
        if any(isnan(temp1))
		  temp = badval;
          temp1(:) = badvalraw;
          temp2(:) = badvalraw;
        end
		%stimco = stimcounter(zzz,ccc,mod2(p,2),a1.stimclassrec(ttt)) + 1;
        data(:,zzz,ccc,mod2(p,2), ...
             stimcounter(zzz,ccc,mod2(p,2),a1.stimclassrec(ttt)) + 1, ...   % a1.stimclassrec tells us the stim number (1-24)
             a1.stimclassrec(ttt)) = temp2;
        [psdvar,f] = pwelch(temp1,hamming(1000),0,2 ^ nextpow2(fsorig),fsorig);
        psd(zzz,ccc,mod2(p,2), ...
             stimcounter(zzz,ccc,mod2(p,2),a1.stimclassrec(ttt)) + 1, ...   % calculating psd from raw data
             a1.stimclassrec(ttt),:) = psdvar(1:fupper,1)';
		stimcounter(zzz,ccc,mod2(p,2),a1.stimclassrec(ttt)) = ...
          stimcounter(zzz,ccc,mod2(p,2),a1.stimclassrec(ttt)) + 1;
      end
        %figureprep([100 100 900 300],1);plot(f,10*log10(psdvar));
        %figurewrite(sprintf('~/psd/task%d_%03d%02d%d', ...
        %    mod2(p,2),ccc,a1.stimclassrec(ttt),stimco));
       
    end
  end
  
  
end


%% CHANNEL AVG TIME AND FREQUENCY RESPONSE

bb_data_avg = mean(reshape(bb_data,length(epochtime_bb),numchannels,numtasks*numreps*numstimuli),3)';
data_avg = mean(reshape(data,length(epochtime),numchannels,numtasks*numreps*numstimuli),3)';
psd_avg = squeeze(mean(reshape(psd,numchannels,numtasks*numreps*numstimuli,fupper),2));


%% PLOT CHANNEL AVG TIME AND FREQUENCY RESPONSE

figure,imagesc(bb_data_avg);
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

figure,imagesc(psd_avg);
colorbar;



% under review below this %


%% CAR
tempsum = 0;
for ccc = 1:96
    tempsum = tempsum + recdata{1,ccc};
end
tempsum = tempsum/96;
car_recdata = recdata;
for ccc = 1:96
    car_recdata{1,ccc} = car_recdata{1,ccc} - tempsum;
end
    
%% AVG ACROSS TRIALS FOR A CHANNEL
psdmean_trials= zeros(length(subjects),numchannels,numtasks,numstimuli-1,fupper,'single');
datamean_trials = zeros(length(epochtime),length(subjects),numchannels,numtasks,numstimuli,'single');


% Data
for lll = 1:numstimuli
    tempsum = 0;
    for rrr = 1:numreps
        tempsum = tempsum + squeeze(data(:,:,:,:,rrr,lll));
    end
    datamean_trials(:,:,:,:,lll) = tempsum/numreps;
end

% PSD
for lll = 1:numstimuli
    tempsum = 0;
    for rrr = 1:numreps
        tempsum = tempsum + squeeze(psd(:,:,:,rrr,lll,:));
    end
    psdmean_trials(:,:,:,lll,:) = tempsum/numreps;
end

% Normalizing PSDs
psdmean_norm = psdmean_trials;

for p = 1:numtasks
    tempsum = 0;
    for lll =1:numstimuli
        tempsum = tempsum + squeeze(psdmean_norm(:,:,p,lll,:));
    end
    tempsum = tempsum/22;
    for lll = 1:numstimuli
        psdmean_norm(:,:,p,lll,:) = log(squeeze(psdmean_norm(:,:,p,lll,:))./(tempsum));
    end
end

%% PCA
total_psd = zeros(numtasks,numchannels*numstimuli,fupper);

counter = 1;
for p = 1:numtasks
    for ccc = 1:numchannels
        for lll = 1:numstimuli
            for ppp = 1:fupper
                total_psd(numtasks,counter,ppp) = [squeeze(psdmean_norm(1,ccc,p,lll,ppp))]';
                if(isinf(total_psd(numtasks,counter,ppp)) | isnan(total_psd(numtasks,counter,ppp)))
                    total_psd(numtasks,counter,ppp)=0;
                end
            end
            counter = counter + 1;
        end
    end
end

[coeff,score,latent,~,explained,mu] = pca(squeeze(total_psd(2,:,:))');

%% PCA plots
score_length = size(score, 1);
tx_vect= score(1:score_length, 1);
plot(tx_vect,'r'); hold on;
tx_vect= score(1:score_length, 2);
plot(tx_vect,'b');
tx_vect= score(1:score_length, 3);
plot(tx_vect,'g'); 
xlim([0 200]); ylim([-10 10]); hold off;


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

figureprep([100 100 900 300],1); hold on; plot(log(smooth(bbtemp{1,46},2000)));straightline(theseOnsets,'v','m-');straightline(theseOnsets+0.1*fsorig,'v','g-');  %subject YBA
% for ccc = 1:numchannels
%     figureprep([100 100 900 300],1); hold on; plot(recdata{1,ccc}); straightline(onsetpoints,'v','m-');hold off;
%     figurewrite(sprintf('~/inout/rawdata/channel%03d',ccc));
% end
figureprep([100 100 900 300],1); plot(bb_bptemp{1,46});
figureprep([100 100 900 300],1); hold on; plot(recdata{1,46});
xlim([15000 30000]);straightline(theseOnsets,'v','m-');straightline(theseOnsets+0.1*fsorig,'v','g-');
figureprep([100 100 900 300],1);plot(10*log10(squeeze(psdmean_trials(1,46,2,4,:))));xlim([0 200])
figurewrite('bb_subjYBA_ch46');
figureprep([100 100 900 300],1); plot(log(smooth(bbtemp{1,78},2000)));  %subject YBA
figurewrite('bb_subjYBA_ch78');
%figureprep([100 100 900 300],1); plot(log(smooth(bbtemp{3,86},2000)));
%figureprep([100 100 900 300],1); plot(log(smooth(bbtemp{3,108},2000)));

figure;plot(mean(catcell(2,bbtemp(1,:)),2));
figurewrite('mean');
%figure;plot(mean(catcell(2,bbtemp(3,:)),2));


for ccc = 1:numchannels
    semilogy((squeeze(psdmean_trials(1,ccc,2,4,:))),'b');
    hold on;
end

%channel maps
for p = 1:numtasks
    for lll = 1:numstimuli
        figureprep([100 100 900 300],1);
		clims = [-250 250];
        imagesc((squeeze(datamean_trials(:,1,:,p,lll)))',clims);
        colorbar;
        set(gca,'XTick',0:1000:8000)
        set(gca,'XTickLabel',-0.5:0.5:3.5)
        figurewrite(sprintf('task%d_stim%02d',p,lll),[],[],'~/inout/channelmaps');
    end
end

%RTPO_3_4 plots
for p = 1:numtasks
    for lll = 1:numstimuli
        figureprep([100 100 1800 300],1);
        subplot(2,2,[1,2]);
        plot((squeeze(datamean_trials(:,1,73,p,lll)))','r');
        hold on;
        plot((squeeze(datamean_trials(:,1,74,p,lll)))','g');
        hold off;
        set(gca,'XTick',0:1000:8000)
        set(gca,'XTickLabel',-0.5:0.5:3.5)
        subplot(2,2,3);
        semilogy((squeeze(psdmean_trials(1,73,p,lll,:))),'r');
        hold on;
        semilogy((squeeze(psdmean_trials(1,74,p,lll,:))),'g');
        hold off;
        subplot(2,2,4);
        plot((squeeze(psdmean_norm(1,73,p,lll,:))),'r');
        hold on;
        plot((squeeze(psdmean_norm(1,74,p,lll,:))),'g');
        hold off;
        figurewrite(sprintf('task%d_stim%02d',p,lll),[],[],'~/inout/RTPO_3_4');
    end
end

% psdmean_trials Plot
semilogy((squeeze(psdmean_trials(1,76,2,3,:))),'m');
hold on;
semilogy((squeeze(psdmean_trials(1,76,2,18,:))),'g');
semilogy((squeeze(psdmean_trials(1,76,2,19,:))),'b');
semilogy((squeeze(psdmean_trials(1,76,2,5,:))),'r');
hold off;
lgd = legend('3','5','8','100');
title(lgd,'Word Contrast')

% psdmean_norm plot
figureprep([100 100 900 300],1);
hold on;
plot((squeeze(psdmean_norm(1,76,2,17,:))),'m');
plot((squeeze(psdmean_norm(1,76,2,18,:))),'g');
plot((squeeze(psdmean_norm(1,76,2,19,:))),'b');
plot((squeeze(psdmean_norm(1,76,2,5,:))),'r');
hold off;
lgd = legend('3','5','8','100');
title(lgd,'Word Contrast')

% OLD:
%         % visualize for sanity
%         figureprep([100 100 1000 300]); hold on;
%         plot(bb);
%         figurewrite(sprintf('ch%d',chantoload),[],[],sprintf('~/inout/bb_subj%d/file%d',zzz,p));

%%
tempsum = recdata{1,47};
for i = 1: size(tempsum,1)
    if(isinf(tempsum(i,1)) | isnan(tempsum(i,1)))
                    tempsum(i,1)=0;
    end
end
recdata{1,47} = tempsum;