%% Power Spectral components (Fig. S2c,d) 
% 
% Performing PCA on PSD to identify freuency bands that explain most variance in the data
%
% ZQ 2023 Mayo Clinic
%
%% PCA analysis

% We compute psd separately in 4 temporal windows (0:0.5:2 s)
% This is an independent script so setup and load data

clc; clear; close all;

dataPath = setLocalDataPath(1);
load(fullfile( dataPath.input, 'derivatives/psd_timechunks.mat'))

bd_tr = [ 30,  54, 136, ... 
         222, 246, 328, ...
         414, 438, 520, ...
         606, 630, 712];

% temp_dir = fullfile( dataPath.output, 'eps/psc');
% if not(isfolder(temp_dir))
%     mkdir(temp_dir)
% end

%% for 3 example electroedes

elec_consider{1} = [75];     % EVC: 75 , VTC: 46, left: 11
elec_consider{2} = [110];    % LVC: 110, VTC: 86, left: 20

for zzz = 1:2

    for ccc = elec_consider{zzz}

        fprintf('\nSub:%d, Ch:%03d',zzz,ccc)

        test_psd = [reshape( psd_on( zzz, ccc, :, :, 4:19, 6:200, 1), [], 195); 
                    reshape( psd_on( zzz, ccc, :, :, 4:19, 6:200, 2), [], 195); 
                    reshape( psd_on( zzz, ccc, :, :, 4:19, 6:200, 3), [], 195);
                    reshape( psd_off(zzz, ccc, :, :, 4:19, 6:200)   , [], 195)];

        test_psd(bd_tr, :) = nan;

        % mean-center and normalize
        avg_psd = mean(test_psd, 1, 'omitnan');
        psd_norm = log10( test_psd) - log10( avg_psd);

        % PCA

        [coeff,score,latent,tsquared,explained] = pca(psd_norm, 'Centered', false);

        % Remove bands with notch (60 harmonics)
        coeff(48:58  , :) = nan;
        coeff(108:121, :) = nan;
        coeff(173:182, :) = nan;

        % Plot frist 2 PCs and their explained variance

        figure,
        hold on
        for ii = 1:2
            plot([6:200], coeff(:,ii), "LineWidth", 2)
        end
        yline(0, 'k--')

        ylim([-0.15 0.25])

        xlabel('Frequency(Hz)')
        ylabel('log-Power')
        title(sprintf('Sub:%d, Ch:%s',zzz,channellabels{zzz}{ccc}))

        legend(sprintf('First PC (%0.2f)', explained(1)), sprintf('Second PC (%0.2f)', explained(2)))%  , sprintf('Third PC (%0.2f)', explained(3)))
        %exportgraphics(gca,[temp_dir sprintf('/Sub%d_el%03d.eps', zzz, elec_consider{zzz})], 'BackgroundColor','none','ContentType','vector', 'Resolution',300);

    end

end