%% Electrode rendering (Fig. 3a,b) 
% 
% Plot electrodes in fsAverage pial and show average HFB power for visual-responsive electrodes
%
% ZQ 2026 Mayo Clinic
%
%% Normalized Power for all visual-reponsive electrodes

stimgroups  = [4, 5];
stimgrnames = {'Word' 'Face'};

% Time window for estimating responsiveness 0 - 400 ms (maches the window length of baseline)
winlen_resp = find(epochtime_bb > 0 & epochtime_bb < 0.4); 

elec_weights = {};

% elec_ROI: visual-responsive electrodes from preprocessing_2_respElec.m

for zzz = 1:2

    % geomean across all the [tasks x trials]

    tempx = reshape( hfbdata_pc(winlen_resp, zzz, elec_ROI{zzz}, :, :, stimgroups), length(winlen_resp), length(elec_ROI{zzz}), []);
    tempx = geomean( geomean( tempx + 1 , 3, 'omitnan'), 1, 'omitnan') - 1;
    elec_weights{zzz} = tempx;

end

clear tempx

% Normalize wrt max fro visualizaiton

temp_mx = ceil( max( [elec_weights{1},  elec_weights{2}]));

tempx = [elec_weights{1},  elec_weights{2}] ./ temp_mx;
assert( length(tempx) ==  length(elec_ROI{1}) + length(elec_ROI{2}))

elec_weights{1} = tempx( 1:length(elec_ROI{1}) );
elec_weights{2} = tempx( length(elec_ROI{1})+1 : end);

clear temp*

%% Run below for both the subjects

addpath(fullfile('/Users/xxx/Documents/git_repo/spm12'));
addpath(fullfile('/Users/xxx/Documents/git_repo/fieldtrip/external/freesurfer/'));

% sub-ID 1 or 2
temp_sub = 2;

electrodes_tsv_name = fullfile( dataPath.input, 'ecogTaskDemands-raw', subjects{temp_sub}, 'ses-ieeg01/iEEG', [subjects{temp_sub} '_ses-ieeg01_space-MNI305_electrodes.tsv']);

loc_info = readtable(electrodes_tsv_name,'FileType','text','Delimiter','\t','TreatAsEmpty',{'N/A','n/a'});

% extract electrode name, xyz co-ordinates, and l/r hemi
% Using std (fsAverage/MNI305) for comparision across 2 subjects
elname = [loc_info.name];
xyz_fsavg = [loc_info.x loc_info.y loc_info.z];
el_hemi   = [loc_info.hemisphere];

hemi = 'R';

views_plot = [-15 -30];  % or [15 30] for left  [-15 30] for left

figure
% load gifti and sulcal labels
g = gifti( fullfile( dataPath.fs, 'fsaverage', ['pial.' lower(hemi) '.surf.gii']));
sulcal_labels = read_curv( fullfile( dataPath.fs, 'fsaverage', 'surf', [lower(hemi) 'h.sulc']));

% plot pial
tH = ieeg_RenderGifti( g, sulcal_labels);
set(gcf, 'CurrentAxes', tH.Parent)

% Find electrodes in this hemi
electrodes_thisHemi = find(ismember(loc_info.hemisphere, hemi));

els_pop = xyz_fsavg(electrodes_thisHemi,:);

% Find visual-responsive electrodes
idx1 = ismember(electrodes_thisHemi,elec_ROI{temp_sub});
idx2 = ismember(elec_ROI{temp_sub},electrodes_thisHemi);
ieeg_elAdd_sizable( els_pop(idx1,:), (elec_weights{temp_sub}(idx2))', 1, 70)
ieeg_elAdd( els_pop(~idx1,:), [.99 .99 .99], 40)    % add electrode positions - white circle on black dot
ieeg_elAdd( els_pop(~idx1,:), 'k', 30)              % add electrode positions - black dot

clear idx1 idx2

% show el names at els
% for ii = 1: length(electrodes_thisHemi)
%     text(xyz(electrodes_thisHemi(ii),1), xyz(electrodes_thisHemi(ii),2), xyz(electrodes_thisHemi(ii),3), elname{electrodes_thisHemi(ii)});
% end

ieeg_viewLight(views_plot(1), views_plot(2)) % change viewing angle
tH.FaceAlpha = 0.9;
set(gcf,'PaperPositionMode','auto')

rmpath(fullfile('/Users/xxx/Documents/git_repo/spm12'));
rmpath(fullfile('/Users/xxx/Documents/git_repo/fieldtrip/external/freesurfer/'));


