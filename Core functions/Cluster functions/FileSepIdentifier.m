function [RootPath,FS]=FileSepIdentifier(RunOnCluster)
% identify the filesep adn main path for cluster or PC

if ismac
    RootPath='/Volumes/';
    FS='/';
elseif ispc
    RootPath='Z:\Projects\Rule_Representation\ElecPhys_Analysis\Rule Representation Project\Submission Code\Submission to Nature March 2024\Tafazoli et al 2024 Code\';
    FS='\';
end
