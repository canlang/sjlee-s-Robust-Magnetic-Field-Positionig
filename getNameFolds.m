function nameFolds = getNameFolds(pathFolder)
        d = dir(pathFolder);
        isub = [d(:).isdir]; %# returns logical vector
        nameFolds = {d(isub).name}';
        nameFolds(ismember(nameFolds,{'.','..'})) = [];
        % nameFolds(~cellfun(@isempty,regexp(nameFolds, '\d{6}'))) = [];
end