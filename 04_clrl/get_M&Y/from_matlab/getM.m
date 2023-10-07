function M = getM(file_path,ubm_path)
    % filePaths = getfilepath(file_path);
    filePaths = file_path;
    num = size(filePaths,1);
    ubm_model = importdata(ubm_path);
    M = [];
    for i = 1:num
        voice_path = filePaths{i,1}
        Mi = spkEnroll(voice_path , ubm_model);
        M = [M,Mi];
    end
    
