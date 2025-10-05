function newNum = mapRef(oldNum)
    % oldNum: scalar or vector of old reference numbers
    % newNum: corresponding new reference numbers (NaN if removed)

    % References that were removed
    removed = [8, 17, 18, 34, 35, 50, 54, 61, 66];

    % Make lookup table (old -> new)
    oldList = 1:68;
    newList = NaN(size(oldList)); % start with NaN for removed
    
    % Fill in new numbering
    counter = 0;
    for i = 1:68
        if ~ismember(i, removed)
            counter = counter + 1;
            newList(i) = counter;
        end
    end

    % Map requested input
    newNum = newList(oldNum);
end