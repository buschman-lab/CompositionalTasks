function opts = ParseOptionalInputs(opts,InputArgs)

%Process optional inputs
if ~isempty(InputArgs)
    if mod(length(InputArgs), 2) ~= 0, error('Must pass key/value pairs for options.'); end
    for i = 1:2:length(InputArgs),
        try
            opts.(InputArgs{i}) = InputArgs{i+1};
        catch
            error('Couldn''t set option ''%s''.', InputArgs{2*i-1});
        end
    end


end

