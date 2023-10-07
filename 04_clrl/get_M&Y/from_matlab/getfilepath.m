function filePaths = getfilepath(folderPath)
    path = fullfile(folderPath);
    files = dir(fullfile(path,'*.wav'));
    fileNames = {files.name}';
    lengthNames = size(fileNames,1);
    filePaths = [];
    for i = 1:lengthNames
        filePath = strcat(path,fileNames(i));
        filePaths = [filePaths;filePath];
    end
