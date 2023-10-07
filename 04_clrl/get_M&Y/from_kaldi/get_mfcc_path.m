function filesPath = get_mfcc_path(folderPath)
    % ��ȡ�ļ���������.mfcc�ļ�
    path = fullfile(folderPath);
    files = dir(fullfile(path,'*.mfcc'));
    filesName = {files.name}';  % ��ȡpath·�������к�׺Ϊmfcc���ļ���
    filesNum = size(filesName,1);  % .mfcc�ļ���
    filesPath = zeros(filesNum,1); % �����ļ�·��
    filesPath = string(filesPath);  % תΪ�ַ�������
    for i = 1:filesNum
        filePath = strcat(path,filesName(i));
        filesPath(i,:) = filePath;
    end
end

