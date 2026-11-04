%% Preprocessing_2 
% 
% Find visual-responsive channels and group into EVC, VTC - face/word selective
%
% ZQ 2023 Mayo Clinic
%
%% Finding responsive channels

% only considering 100% contrast word/face stim
stimgroups  = [4,5];
stimgrnames = {'Word' 'Face'};

p = NaN(length(subjects), numchannels);
q = NaN(length(subjects), numchannels);
respChs = {};    % responsive channels

% Baseline time window -500 to -100 ms
winlen_base = find(epochtime_bb >= -0.5 & epochtime_bb < -0.1);

% Time window for estimating responsiveness 0 - 400 ms (maches the window length of baseline)
winlen_resp = find(epochtime_bb > 0 & epochtime_bb <= 0.4); 


for zzz = 1:2 % subjects
    for ccc = gcc{zzz}  % recording channels

        % combining stim x trials x task, then geomean

        temp_resp = reshape( hfbdata_pc( winlen_resp, zzz, ccc, :, :, stimgroups), length(winlen_resp), []);
        x = geomean( temp_resp+1, 1, 'omitnan') - 1;
        
        temp_base = reshape( hfbdata_pc( winlen_base, zzz, ccc, :, :, stimgroups), length(winlen_base), []);
        x_base = geomean( temp_base+1, 1, 'omitnan') - 1;
        
        % Wilcoxon signed rank test - test for median
        p (zzz, ccc) = signrank(x, x_base, 'tail', 'right', 'Alpha', 0.05);

    end 

    % use BHFDR corrrection for multiple comparision
    val_idx = ~isnan(p(zzz, :));
    q(zzz, val_idx) = mafdr(p(zzz, val_idx), 'BHFDR', true);

    % visual-responsive channels at the 5% significance level BHFDR corrected
    respChs{zzz} = find( q(zzz,:) <= 0.05 );

end


clear x x_base temp* val_idx zzz ccc 


%% Manually select channels in ROIs of interest: EVC & VTC

elec_ROI = {[9, 10, 11, 12, 15, 41, 44, 45, 46, 72, 74, 75, 76, 77, 78], ...    % sub1
            [17, 20, 21, 22, 82, 84, 86, 87, 88, 106, 108, 110]};           % sub2  

elec_EVC = {[72, 74, 75, 76], [106, 108]};
elec_VTC = {[9, 10, 11, 12, 15, 41, 44, 45, 46, 77, 78], [17, 20, 21, 22, 82, 84, 86, 87, 88, 110]};

%% Category Selectivity using d-prime
% categorizing VTC channels into face-, word-, and non-selective - VTC based on d' values

dprime = NaN(length(subjects), numchannels);
 
for zzz = 1:2 % subjects
    for ccc = elec_VTC{zzz}

        % combining stim x trials x task, then geomean

        temp_resp = squeeze( reshape( hfbdata_pc( winlen_resp, zzz, ccc, :, :, stimgroups), length(winlen_resp), [], 2));
        x = squeeze( geomean( temp_resp+1, 1, 'omitnan'));
        
        % compute d-prime
        temp_meanx = geomean( x, 1, 'omitnan') - 1;
        temp_varx  = exp( var ( log(x), 1, 'omitnan')) - 1;
        dprime (zzz, ccc) = (temp_meanx(1) - temp_meanx(2)) / ...
                            sqrt(0.5 * (temp_varx(1) + temp_varx(2)));     % word - faces
        
    end
end

% sort into categories based on d-prime
for zzz = 1:2
    temp_ws = [];
    temp_fs = [];
    temp_ns = [];
    for ccc = elec_VTC{zzz}
        if dprime (zzz, ccc) >= 1       % word-selective       
            temp_ws = [temp_ws, ccc];
        elseif dprime (zzz, ccc) <= - 1 % face-selective 
            temp_fs = [temp_fs, ccc]; 
        else
            temp_ns = [temp_ns, ccc];   %  non-selective
        end
    end
    elec_VTC_w {zzz} = temp_ws;   % word selective 
    elec_VTC_f {zzz} = temp_fs;   % face selective
    elec_VTC_ns{zzz} = temp_ns;   %  non-selective
end

clear x temp*