function [FileName,Path,FullPath] = GenerateFileName(FS,Path,AnalysisName,Animal,Date,Ch,ExtraParams,varargin)
%GENERATEFILENAME
% generates a filename based on below template
% toggle SelfName to remove the template Don't include '.mat' in it
global AnalysisOpts
%AnalysisName_Animal_Date_Channel_ExtraParams
opts.SelfName=0;
opts.SelfNameTxt='Untitled';
opts.ext=''; % extention of the file
opts.MakeFolder=0; % make a folder with name of the file in the same directory for everything related goes there
opts.ExtFolderName='';
opts.UseReRef=1; % use reref to change the name
% depending on the analysis are we saving the channel names or
% just pulttign all of them toghether
if contains(AnalysisName,'classifier','IgnoreCase',1) | contains(AnalysisName,'subspace','IgnoreCase',1)
    AnalysisOpts.UseIndividualChName=0;
else 
    AnalysisOpts.UseIndividualChName=1;
end
if mod(length(varargin), 2) ~= 0, error('Must pass key/value pairs for options.'); end
for i = 1:2:length(varargin)
    try
        opts.(varargin{i}) = varargin{i+1};
    catch
        error('Couldn''t set option ''%s''.', varargin{2*i-1});
    end
end
% change/adjust aninal name

if ~isempty(AnalysisOpts.DateNum) & (sum(AnalysisOpts.DateNum)==0 | isempty(setdiff(AnalysisOpts.DateSet_2look,AnalysisOpts.DateNum))) % if we are using all of the recordings
    Animal='ALL';Date='ALL';
end
if length(AnalysisOpts.DateNum)>1 
    Date='ALL';
end
if ~isempty(AnalysisOpts.Animal); Animal=AnalysisOpts.Animal;end% update animal 
if ~strcmp(Date,'ALL') & isempty(Animal)
    if strcmpi(AnalysisOpts.Project,'Rule Representation')
        Animal=AnalysisOpts.AnimalSet{strcmp(Date,AnalysisOpts.DateSet)};
    elseif strcmpi(AnalysisOpts.Project,'Learning attentional templates')
        Animal=AnalysisOpts.AnimalSet{strcmp(Date(3:end),AnalysisOpts.DateSet)};
    end    
elseif strcmp(Date,'ALL') & isempty(Animal)
    Animal='ALL';
    Date='ALL';
elseif strcmp(Date,'ALL') & ~isempty(Animal)
    Date='ALL';
end
if strcmp(Date,'ALL') % if Date is ALL then folder is ALL 
    opts.ExtFolderName='ALL';
end
if ~opts.SelfName
    if length(Ch)==1 & AnalysisOpts.UseIndividualChName% it is only channel number
        if   AnalysisOpts.ReRefChannels & opts.UseReRef % add rereference or not
            FileNameNoExt=[AnalysisName '_' Animal '_' Date '_' num2str(Ch) ExtraParams '_' AnalysisOpts.ExtraStr];
%          elseif strcmp(Date,'ALL') & ~AnalysisOpts.UseIndividualChName
%              FileNameNoExt=[AnalysisName '_' Animal '_' Date ExtraParams];       
        else
            FileNameNoExt=[AnalysisName '_' Animal '_' Date '_' num2str(Ch) ExtraParams];
        end
    elseif length(Ch)==2  & AnalysisOpts.UseIndividualChName % then we are adding the cluster number to this channel as well 
        if isnan(Ch(2));Clusttxt='';else;Clusttxt=['_cl' num2str(Ch(2))];end
        if   AnalysisOpts.ReRefChannels & opts.UseReRef % add rereference or not
            FileNameNoExt=[AnalysisName '_' Animal '_' Date '_' num2str(Ch(1)) Clusttxt ExtraParams '_' AnalysisOpts.ExtraStr];
%          elseif strcmp(Date,'ALL') & ~AnalysisOpts.UseIndividualChName
%              FileNameNoExt=[AnalysisName '_' Animal '_' Date ExtraParams];       
        else
            FileNameNoExt=[AnalysisName '_' Animal '_' Date '_' num2str(Ch(1)) Clusttxt ExtraParams];
        end
    elseif strcmp(Ch,'ALL') |  ~AnalysisOpts.UseIndividualChName
        FileNameNoExt=[AnalysisName '_' Animal '_' Date ExtraParams];
    end
    FileNameNoExt=strrep(FileNameNoExt,'.',''); % remove any dot from string
    FileName=[FileNameNoExt AnalysisOpts.ExtaStrAnalysis opts.ext]; % add extra string if we need to at last moment 
else
    [~,FileNameNoExt]=fileparts(opts.SelfNameTxt);
    FileName=[opts.SelfNameTxt];
end
if AnalysisOpts.RunDummyFile;FileName=['Dum_' FileName];end % if we are running a dummy file 
    
if opts.MakeFolder % make a folder for this file
    % Path=[Path AnalysisName FS FileNameNoExt FS];
    Path=[Path AnalysisName FS opts.ExtFolderName FS]; % create a subfolder with name of this analysis
else
    if isempty(Date) & isempty(opts.ExtFolderName)
        Path=[Path AnalysisName FS];
    elseif ~isempty(Date) & isempty(opts.ExtFolderName)
        Path=[Path AnalysisName FS Date FS];
    elseif ~isempty(opts.ExtFolderName)
        Path=[Path AnalysisName FS opts.ExtFolderName FS];
    end
end

if ~exist(Path,'file')
    mkdir(Path)
end

FullPath=[Path FileName];

end

