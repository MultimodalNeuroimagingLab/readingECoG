%% Spectral power analysis (Fig. 3c,d) 
% 
% Spectrograms for visual-responsive electrodes
%
% ZQ 2023 Mayo Clinic
%
%% Example electrodes to plot

elec_consider{1} = [75];      % EVC: 75 , VTC: 46, left: 11
elec_consider{2} = [110];     % EVC: 110, VTC: 86, left: 20


% temp_dir = fullfile( dataPath.output, 'eps/spectra');
% if not(isfolder(temp_dir))
%     mkdir(temp_dir)
% end 

%% Spectograms using Multi-tapper 

% add chronux_2_11 to the path
% Note this might interfere with other matlab in-built functions like findpeaks
addpath(genpath('/Users/xxx/Documents/Matlab/chronux_2_11/spectral_analysis'))

% Define time-frequency muti-tapper settings
movingwin=[.300 .01];
params.pad=-1;
params.tapers=[5 5];
params.fpass=[0 200];
params.Fs=fsorig;
params.trialave=0;

% baseline window: -500 to -100 ms
winlen_base = find(epochtime >= -0.5 & epochtime < -0.1);
spect_base = nan(length(subjects),numchannels,260);

stimgroups  = {5, 4};    % only Face 100% Con/PC response for sub 1 and Word 100% Con/PC response for sub 2
stimleg  = {'face', 'word'};

% response window fro visualization
winlen_tc = find(epochtime >= -0.5 & epochtime <= 2.5);   % -0.5 s to 2.5 s

% subject 1 & 2
for zzz = 1:2
    for ccc =  elec_consider{zzz}

       % merge across task x trials x stimulus for baseline
        temp_base = reshape( data( winlen_base, zzz, ccc, :, :, :), length(winlen_base), []);

        temp_base = rmmissing(temp_base')';     % remove NaNs
        
        % calculate specgram
        [temp_spect, ~, ~] = mtspecgramc(temp_base, movingwin, params);

        % median across trials and time
        spect_base = median( median( temp_spect, 3, 'omitnan'),  1, 'omitnan');


        % merge across task x trials for responses
        temp_signal = reshape( data( winlen_tc, zzz, ccc, :, :, stimgroups{zzz}), length(winlen_tc),[]);
        
        % calculate specgram
        [temp_spect, t_spect, f_spect]= mtspecgramc(temp_signal, movingwin, params);
        
        % median across trials
        temp_spect = median( temp_spect, 3, 'omitnan');
        
        % make sure the first time point is correct
        t_spect = t_spect + epochtime(1);
        
        % log-normalize wrt baseline
        spectra = log10(temp_spect) - log10(spect_base);
        
        figure
        load('external/dg_colormap.mat', 'cm') % KJM color scheme;

        uimagesc(t_spect,f_spect,spectra', [-1 1]); axis xy; colormap(cm); colorbar
        ylabel('Frequency (Hz)');
        xlabel('Time (s)');
        title(sprintf('Sub:%d, Ch:%s',zzz,channellabels{zzz}{ccc}))

        %exportgraphics(gcf,[temp_dir sprintf('/Subj%d_el%d_%s.eps',zzz,ccc, stimleg{zzz})], Resolution=300);

    end
end

clear temp*

rmpath(genpath('/Users/xxx/Documents/Matlab/chronux_2_11/spectral_analysis'))
