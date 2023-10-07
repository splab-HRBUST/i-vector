function filesPath = get_mfcc_path(folderPath)
    % 读取文件夹下所有.mfcc文件
    path = fullfile(folderPath);
    files = dir(fullfile(path,'*.mfcc'));
    filesName = {files.name}';  % 获取path路径下所有后缀为mfcc的文件名
    filesNum = size(filesName,1);  % .mfcc文件数
    filesPath = zeros(filesNum,1); % 保存文件路径
    filesPath = string(filesPath);  % 转为字符串类型
    for i = 1:filesNum
        filePath = strcat(path,filesName(i));
        filesPath(i,:) = filePath;
    end
end

