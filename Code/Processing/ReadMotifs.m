function [motif_onset,st_norm] = ReadMotifs(cur_rec,motif)
if ispc
    addpath(genpath('Z:\Rodent Data\Wide Field Microscopy\fpCNMF'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Widefield_Imaging_Analysis'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\GenUtils'));
    addpath(genpath('Z:\Projects\Cortical Dynamics\Cortical Neuropixel Widefield Dynamics\GithubRepo\Ephys'));
    savedir = 'Z:\Projects\Oscillation in Cortical Dynamics\Input Output Data\';
else
    addpath(genpath('/jukebox/buschman/Rodent Data/Wide Field Microscopy/fpCNMF'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Widefield_Imaging_Analysis'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/GenUtils'));
    addpath(genpath('/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/GithubRepo/Ephys'));
    savedir = '/jukebox/buschman/Projects/Cortical Dynamics/Cortical Neuropixel Widefield Dynamics/analysisplayground/CCA/';
end
tic
win=[-2 9]; %hardcoded write now.
%starting
fprintf('Working on rec %d motif %d',cur_rec,motif);
%% Gathering Data
[rec_name,~,~,EphysPath,motif_fits] = LoadDataDirectories(cur_rec);
%load ephys data
[st_mat,~,st_depth] = LoadSpikes(EphysPath,'bindata',1,'offset',15,'mua',1,'depth_type','probe');
% st_norm = cellfun(@(x) x./std(x),st_mat,'UniformOutput',0);
st_norm = st_mat;
st_norm = cellfun(@(x) x(1:2:end,:)+x(2:2:end,:),st_norm,'UniformOutput',0);
%get the anatomical locations
neu_area = LoadEphysAnatomicalLocations(EphysPath,st_depth);
neu_area = cat(2,neu_area{:});
%get the motif onsets for all recordings
[motif_onset,~] = CompileMotifOnsets(motif_fits); %return 'chunk_names' if you want to confirm ordering is correct
%motif_onset = cellfun(@(x) floor(x/2), motif_onset,'UniformOutput',0);

% if motif > numel(motif_onset) %get null periods that of window length
%     fprintf('\n\t Running on NULL')
%     [~,trig_st] = ParseNullPeriods([],st_norm,motif_onset,win,10);
% else %parse motif onsets
%     fprintf('\n\t Running on Motifs')
%     [~,trig_st] = ParseByOnset([],st_norm,motif_onset,win,motif);
% end
% %parse activity per parent region
% [area_val, area_label] = ParseByArea(cat(1,trig_st{:}),neu_area,'general');

end

