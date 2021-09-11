clear all
dataRootPath = uigetdir;    %'//Users/zeeshanqadir/Documents/git_repo/readingECoG/readingECoG-raw/sub-YBA/ses-ieeg01/ieeg/sub-YBA_ses-ieeg01_task-categorization_run-01_ieeg';

prompt = {'Enter subject label:','Enter session label:','Enter task label:','Enter run label:'};
dlgtitle = 'Input';
answer = inputdlg(prompt,dlgtitle);

sub_label = answer{1};  %'01';
ses_label = answer{2};  %'01';
task_label = answer{3};  %'fixation';
run_label = answer{4};  %'01';
 
ieeg_name = fullfile(dataRootPath,['sub-' sub_label],['ses-' ses_label],'ieeg',...
    ['sub-' sub_label...
    '_ses-' ses_label...
    '_task-' task_label...
    '_run-' run_label ...
    '_ieeg.mat']);
 
ieeg_channels_name = fullfile(dataRootPath,['sub-' sub_label],['ses-' ses_label],'ieeg',...
    ['sub-' sub_label...
    '_ses-' ses_label...
    '_task-' task_label...
    '_run-' run_label ...
    '_channels.tsv']);
 
ieeg_name_save = fullfile(dataRootPath,['sub-' sub_label],['ses-' ses_label],'ieeg',...
    ['sub-' sub_label...
    '_ses-' ses_label...
    '_task-' task_label...
    '_run-' run_label ...
    '_ieeg']);
 
% load the data
disp(['loading raw data ' ieeg_name])
load(ieeg_name);
data = data';
 
% load channel info 
disp(['loading channels.tsv ' ieeg_channels_name])
channel_info = readtable(ieeg_channels_name,'FileType','text','Delimiter','\t');
%%              
 
% data = data(:,1:3000);
 
% make chan_labels a cell array:        
if ~ischar(channel_info.name{1})
    chan_labels = textscan(num2str(channel_info.name'),'%s');
    chan_labels = chan_labels{1};
else
    chan_labels = channel_info.name;
end
        
% assign header fields:
dataStruct.hdr.Fs = round(channel_info.sampling_frequency(1));
dataStruct.hdr.nChans = size(data,1);       
dataStruct.hdr.label = chan_labels;
dataStruct.hdr.nSamples = size(data,2);
dataStruct.hdr.nSamplesPre = 0;
dataStruct.hdr.nTrials = 1;
dataStruct.hdr.chantype = channel_info.type;
dataStruct.hdr.chanunit = repmat({'uV'},size(data,2),1);
        
% assign other fields:
dataStruct.label = chan_labels;
dataStruct.time{1} = [1:dataStruct.hdr.nSamples]/dataStruct.hdr.Fs;
dataStruct.trial{1} = data;
dataStruct.fsample = dataStruct.hdr.Fs; % this needs to be a scalar
dataStruct.sampleinfo = [1 dataStruct.hdr.nSamples];
 
hdr_data = ft_fetch_header(dataStruct);
 
% ft_write_data(ieeg_name_save,dataStruct.trial{1},'header',hdr_data,'dataformat','gdf')
ft_write_data(ieeg_name_save,dataStruct.trial{1},'header',hdr_data,'dataformat','brainvision_eeg')
% ft_write_data(ieeg_name_save,dataStruct.trial{1},'header',hdr_data,'dataformat','edf')
        
%% test and load data
 
test_data      = ft_read_data([ieeg_name_save '.eeg'],'dataformat','brainvision_eeg');
test_header    = ft_read_header([ieeg_name_save '.vhdr'],'dataformat','brainvision_eeg');
 
length(find(test_data-dataStruct.trial{1}~=0))