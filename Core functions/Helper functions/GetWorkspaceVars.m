
% Display and sort variables in the workspace by size in descending order
workspace_vars = whos;
[~, idx] = sort([workspace_vars.bytes], 'descend');
sorted_vars = workspace_vars(idx);

for i = 1:length(sorted_vars)
    var_name = sorted_vars(i).name;
    var_size_bytes = sorted_vars(i).bytes;
    var_size_MB = var_size_bytes / 1e6; % Convert bytes to megabytes
    fprintf('Variable: %s, Size: %.2f MB\n', var_name, var_size_MB);
end
