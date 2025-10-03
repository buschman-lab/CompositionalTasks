% Specify the folder containing the .mat files
folderPath = '/path/to/your/folder';

% Specify the variable to be deleted
variableToDelete = 'yourVariableName';

% Get a list of all .mat files in the folder
matFiles = dir(fullfile(folderPath, '*.mat'));

% Loop through each .mat file
for i = 1:length(matFiles)
    % Load the mat file
    filePath = fullfile(folderPath, matFiles(i).name);
    data = load(filePath);
    
    % Check if the variable exists in the loaded data
    if isfield(data, variableToDelete)
        % Remove the specified variable
        data = rmfield(data, variableToDelete);
        
        % Create a new filename (e.g., appending '_modified' to the original filename)
        [~, fileName, fileExt] = fileparts(matFiles(i).name);
        newFileName = fullfile(folderPath, [fileName, '_modified', fileExt]);
        
        % Save the modified data to the new file
        save(newFileName, '-struct', 'data');
        
        % Delete the original file
        delete(filePath);
        
        fprintf('Variable "%s" deleted from file: %s\n', variableToDelete, matFiles(i).name);
        fprintf('New file created: %s\n', newFileName);
    else
        fprintf('Variable "%s" not found in file: %s\n', variableToDelete, matFiles(i).name);
    end
end
