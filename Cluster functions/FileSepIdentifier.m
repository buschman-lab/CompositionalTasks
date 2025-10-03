function [RootPath,FS]=FileSepIdentifier(RunOnCluster)
% identify the filesep adn main path for cluster or PC

if RunOnCluster
    RootPath='/jukebox/buschman/';
    FS='/';
else
   if ismac
         RootPath='/Volumes/buschman/';
         FS='/';
    elseif ispc
        RootPath='Z:\';
        FS='\';
    end
end