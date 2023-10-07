function M = spkEnroll(path, ubm)
    % 提取iv
    mmm = extract_mfcc(path,'0');  % 提取mfcc -- 60*433
    gmm = mapAdapt(mmm, ubm, 10.0, 'm');  %提取每一个音频文件的GMM模型
    % gmm --> w(1*1024),mu(60*1024),sigma(60*1024)
    M = gmm.mu(:);  % 提取GMM模型的均值超矢量
end
