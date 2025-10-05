function [file_list, folder_list] = GrabFiles(substring, interactiveFlag, searchdir)

% @Synopsis 
% Camden MacDowell 2018 Searches folder or subfolders for any file with with defined substring.
% Alternatively user can select individual files. Outputs a list of file
% paths and folder paths. 

% @Inputs 
% @substring (required)
% A string containing the substring of the files to grab. Will grab any file 
% with substring anywhere in the file name; 
% Example 1: substring = '*.pdf' will grab all pdfs in a folder/subfolders
% Example 2: substring = 'applesauce.mat' will grab all files that
% contain applesauce in the name applesauce irrespective of other 
% text in the file name. 
%
%
% @interactiveFlag (optional). 
% If 1 (def) then user can select folder, subfolder, individual files or
% any combination of the above. If 0 then the function only grab files with
% designated substring in the directories specified in searchdir
%
% @searchdir (optional - only applies if interactiveflag = 0)
% Cell array of directories to search for files if interactiveflag is 0; 
%
% @outputs: list of paths for each file and folder paths for each folder

if nargin < 2
    interactiveFlag = 1; 
end
if nargin <3 
    searchdir = {pwd};     
end

%Grab files 
file_list = [];
folder_list = [];
start = pwd;

if interactiveFlag %if user interaction 
    while 1 %Log selection loop
        dlg = questdlg('Would you like to select all files in a folder, a subfolder, or individual files?', ...
            'File vs Folder vs Subfolder Selection', ...
            'Folders','SubFolders','Files','Folders');
        switch dlg
            case 'Folders' %Select the folder containing logs
                target_dir = uigetdir();
                [file_names,file_path, folder_names] = GrabFile(target_dir,substring); %This grabs just the image stacks, not individual images. 
                file_list = [file_list, file_path];
                folder_list = [folder_list, folder_names];
                cd(start);
            case 'Files' %Pick the log
                [baseName, folder] = uigetfile(substring,'Pick a File');
                file_path = fullfile(folder, baseName);
                folder_list = [folder_list, {folder}];
                file_list = [file_list, {file_path}];
            case 'SubFolders' %Select the folder containing subfolders containing logs
                target_dir = uigetdir();
                tempd = dir(target_dir);
                isub = [tempd(:).isdir]; %returns logical vector
                nameFolds = {tempd(isub).name}';
                nameFolds(ismember(nameFolds,{'.','..'})) = []; %remove silly extra folders
                %Loop through each subfolder in the dir. 
                for i = 1:length(nameFolds) %subfolder loop
                    subfolder_fn = [target_dir filesep num2str(nameFolds{i})];
                    [file_names,file_path,folder_names] = GrabFile(subfolder_fn,substring); 
                    file_list = [file_list, file_path];
                    folder_list = [folder_list,folder_names];
                end %Subfolderloop
                cd(start)

        end %end switch
        choice = questdlg(sprintf('Would you like to select additional %s?',dlg), ...
            'Batch Files', ...
            'Yes','No','Yes');

            % Handle response
            switch choice
                case 'Yes'

                case 'No'
                    break %folder selection while loop
            end
    end %log selection while loop
else %just grab files in the directories listed in searchdir
    for cur_f = 1:size(searchdir)
        target_dir = searchdir{cur_f};
        [file_names,file_path, folder_names] = GrabFile(target_dir,substring); %see subfunc below
        file_list = [file_list, file_path];
        folder_list = [folder_list, folder_names];
        cd(start);
    end
end %interactive if/else

%Remove duplicates if they occur
file_list = unique(file_list);
folder_list = unique(folder_list);
cd(start);


function [in_fns,in_path,in_folder] = GrabFile(target_dir,substring)
if isempty(target_dir) 
    target_dir = uigetdir;
else
end
startdir = (pwd);
cd(target_dir);
files = dir(); 
matches = regexpi({files.name},substring);
files = files(cellfun(@(x)~isempty(x),matches));
% files = dir(substring); %just grab those with a specific substring at the end
if isempty(files)   %No files = throw an error
    warning(sprintf('No files with "%s" in selected directory %s',substring,target_dir));
    in_fns = [];
    in_path = [];
    in_folder = [];
else
    in_fns = cell([1,length(files)]);
    in_path = cell([1,length(files)]);
    in_folder = cell([1,length(files)]);
    for cur_f = 1:length(files)
        in_fns{cur_f} = files(cur_f).name;
        in_folder{cur_f} = files(cur_f).folder;
        in_path{cur_f} = fullfile(files(cur_f).folder,files(cur_f).name);
    end
    in_fns = in_fns(~cellfun('isempty',in_fns));
    in_path = in_path(~cellfun('isempty',in_path));
    in_folder = in_folder(~cellfun('isempty',in_folder));
    if isempty(in_fns)
        warning(sprintf('No files with "%s" in selected directory %s',substring,target_dir));
    end
cd(startdir);
end
