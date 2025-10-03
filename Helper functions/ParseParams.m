function ParseParams(inputArgs)

    global AnalysisOpts
    if ~isempty(inputArgs)
        if iscell(inputArgs{1}) && ~isempty(inputArgs{1}) %length(inputArgs{1})>1 & ~isempty(inputArgs{1}{1})
            if iscell(inputArgs{1}{1})  
                inputArgs=inputArgs{1}{1};
            else
                 inputArgs=inputArgs{1};
            end
        
        elseif isempty(inputArgs{1})
             inputArgs=[];
        end
    end % this is when there is a nested  varargin
    % parse params
    if mod(length(inputArgs), 2) ~= 0, error('Must pass key/value pairs for options.'); end
    for i = 1:2:length(inputArgs)
        try
            indp=findstr(inputArgs{i},'.');
            if length(indp)==1                 
                AnalysisOpts.(inputArgs{i}(1:indp-1)).(inputArgs{i}(indp+1:end))=inputArgs{i+1};
            elseif length(indp)==2
                AnalysisOpts.(inputArgs{i}(1:indp(1)-1)).(inputArgs{i}(indp(1)+1:indp(2)-1)).(inputArgs{i}(indp(2)+1:end))=inputArgs{i+1};                
            else
                AnalysisOpts.(inputArgs{i}) = inputArgs{i+1};
            end
        catch
            error('Couldn''t set option ''%s''.', inputArgs{2*i-1});
        end
    end
end

